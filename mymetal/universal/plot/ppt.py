"""PowerPoint slide export helpers.

Functions:
    ppt2picture: Export every slide in one PowerPoint file to image files.

The conversion uses the Windows PowerPoint COM interface and Pillow. These
are imported inside :func:`ppt2picture` so importing ``mymetal`` remains
possible on systems without Microsoft PowerPoint.
"""

from __future__ import annotations

import gc
from collections import Counter
from pathlib import Path
from typing import Any

from PIL import Image, ImageChops


SUPPORTED_PPT_SUFFIXES = {".ppt", ".pptm", ".pptx"}
SUPPORTED_STYLES = {
    "png": ("PNG", ".png"),
    "jpg": ("JPG", ".jpg"),
    "jpeg": ("JPG", ".jpg"),
    "bmp": ("BMP", ".bmp"),
    "gif": ("GIF", ".gif"),
    "tif": ("TIFF", ".tif"),
    "tiff": ("TIFF", ".tiff"),
}

# PowerPoint files are trusted local inputs and can legitimately export large
# slides. Pillow's default pixel limit would reject some valid figure slides.
Image.MAX_IMAGE_PIXELS = None


def _normalize_style(style: str) -> tuple[str, str, str]:
    """Return the PowerPoint format, output suffix, and Pillow format."""
    if not isinstance(style, str):
        raise TypeError("style must be a string such as 'png' or 'jpg'.")
    style_normalized = style.lower().lstrip(".")
    if style_normalized not in SUPPORTED_STYLES:
        supported = ", ".join(sorted(SUPPORTED_STYLES))
        raise ValueError(
            "Unsupported style: " + style + ". Choose from: " + supported
        )
    powerpoint_format, suffix = SUPPORTED_STYLES[style_normalized]
    pillow_format = "JPEG" if powerpoint_format == "JPG" else powerpoint_format
    return powerpoint_format, suffix, pillow_format


def _check_input(path_ppt: Path) -> Path:
    """Validate and normalize one PowerPoint path."""
    path_ppt = Path(path_ppt).expanduser().resolve()
    if not path_ppt.is_file():
        raise FileNotFoundError("PowerPoint file does not exist: " + str(path_ppt))
    if path_ppt.suffix.lower() not in SUPPORTED_PPT_SUFFIXES:
        raise ValueError(
            "Unsupported PowerPoint suffix: "
            + path_ppt.suffix
            + ". Expected .ppt, .pptm, or .pptx."
        )
    return path_ppt


def _check_parameters(dpi: int, padding_mm: float, crop_threshold: int) -> None:
    """Validate image export parameters."""
    if not isinstance(dpi, int) or isinstance(dpi, bool) or dpi <= 0:
        raise ValueError("dpi must be a positive integer.")
    if padding_mm < 0:
        raise ValueError("padding_mm must be non-negative.")
    if not isinstance(crop_threshold, int) or isinstance(crop_threshold, bool):
        raise ValueError("crop_threshold must be an integer from 1 to 255.")
    if not 1 <= crop_threshold <= 255:
        raise ValueError("crop_threshold must be an integer from 1 to 255.")


def _get_background_color(image: Image.Image) -> tuple[int, int, int]:
    """Estimate the slide background from its four corners."""
    width, height = image.size
    patch = max(4, min(width, height) // 200)
    lpixel: list[tuple[int, int, int]] = []
    lbox = [
        (0, 0, patch, patch),
        (width - patch, 0, width, patch),
        (0, height - patch, patch, height),
        (width - patch, height - patch, width, height),
    ]
    for box in lbox:
        with image.crop(box) as image_patch:
            lpixel.extend(image_patch.getdata())
    return Counter(lpixel).most_common(1)[0][0]


def _crop_image_whitespace(
    path_image: Path,
    dpi: int,
    padding_mm: float,
    crop_threshold: int,
    pillow_format: str,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Crop slide background and save it in the requested image format."""
    with Image.open(path_image) as image_source:
        image_rgb = image_source.convert("RGB")
        original_size = image_rgb.size
        background_color = _get_background_color(image_rgb)
        background = Image.new("RGB", image_rgb.size, background_color)
        difference = ImageChops.difference(image_rgb, background).convert("L")
        mask = difference.point(
            lambda value: 255 if value > crop_threshold else 0
        )
        bbox = mask.getbbox()

        if bbox is None:
            image_result = image_rgb.copy()
        else:
            padding_px = round(dpi * padding_mm / 25.4)
            left = max(0, bbox[0] - padding_px)
            top = max(0, bbox[1] - padding_px)
            right = min(image_rgb.width, bbox[2] + padding_px)
            bottom = min(image_rgb.height, bbox[3] + padding_px)
            image_result = image_rgb.crop((left, top, right, bottom))

    dict_save_kwargs: dict[str, Any] = {
        "format": pillow_format,
        "dpi": (dpi, dpi),
    }
    if pillow_format == "JPEG":
        dict_save_kwargs.update({"quality": 95, "subsampling": 0})
    image_result.save(path_image, **dict_save_kwargs)
    result_size = image_result.size
    image_result.close()
    return original_size, result_size


def _check_exported_image(
    path_image: Path,
    pillow_format: str,
) -> tuple[int, int]:
    """Validate one exported image and return its pixel dimensions."""
    if not path_image.is_file() or path_image.stat().st_size == 0:
        raise RuntimeError("Missing or empty image: " + str(path_image))
    with Image.open(path_image) as image:
        if image.format != pillow_format:
            raise RuntimeError(
                "Unexpected image format for "
                + str(path_image)
                + ": "
                + str(image.format)
            )
        image.verify()
    with Image.open(path_image) as image:
        return image.size


def _get_output_path(
    path_output: Path,
    path_ppt: Path,
    slide_number: int,
    slide_count: int,
    suffix: str,
) -> Path:
    """Build a stable output name for a slide."""
    if slide_count == 1:
        return path_output / (path_ppt.stem + suffix)
    return path_output / (
        path_ppt.stem + "-slide-" + str(slide_number).zfill(2) + suffix
    )


def ppt2picture(
    path_ppt: str | Path,
    path_output: str | Path,
    style: str = "png",
    dpi: int = 600,
    padding_mm: float = 1.0,
    crop_threshold: int = 8,
    overwrite: bool = False,
) -> list[Path]:
    """Export all slides in one PowerPoint file as cropped image files.

    Args:
        path_ppt: Input ``.ppt``, ``.pptm``, or ``.pptx`` file.
        path_output: Directory receiving one image per slide.
        style: Image format, such as ``png``, ``jpg``, ``jpeg``, ``bmp``,
            ``gif``, ``tif``, or ``tiff``. Defaults to ``png``.
        dpi: PowerPoint export resolution in dots per inch. Defaults to 600.
        padding_mm: Background padding retained around slide content, in mm.
        crop_threshold: Pixel difference from the estimated background required
            to count as content, from 1 to 255.
        overwrite: Replace existing output images when True. Existing valid
            images are skipped when False.

    Returns:
        A list of output image paths, including newly exported and skipped files.

    Raises:
        ImportError: If ``pywin32`` is unavailable.
        FileNotFoundError: If the input PowerPoint file is missing.
        RuntimeError: If PowerPoint or image export fails.
    """
    path_ppt = _check_input(Path(path_ppt))
    path_output = Path(path_output).expanduser().resolve()
    path_output.mkdir(parents=True, exist_ok=True)
    powerpoint_format, suffix, pillow_format = _normalize_style(style)
    _check_parameters(dpi, padding_mm, crop_threshold)

    try:
        import pythoncom
        import win32com.client
    except ImportError as error:
        raise ImportError(
            "ppt2picture requires Microsoft PowerPoint and the pywin32 package "
            "on Windows."
        ) from error

    lpath_image: list[Path] = []
    presentation = None
    powerpoint = None
    pythoncom.CoInitialize()
    try:
        powerpoint = win32com.client.DispatchEx("PowerPoint.Application")
        powerpoint.DisplayAlerts = 0
        presentation = powerpoint.Presentations.Open(
            str(path_ppt), True, False, False
        )
        slide_count = int(presentation.Slides.Count)
        width_px = round(float(presentation.PageSetup.SlideWidth) * dpi / 72)
        height_px = round(float(presentation.PageSetup.SlideHeight) * dpi / 72)

        for slide_number in range(1, slide_count + 1):
            path_image = _get_output_path(
                path_output,
                path_ppt,
                slide_number,
                slide_count,
                suffix,
            )
            if path_image.exists() and not overwrite:
                image_size = _check_exported_image(path_image, pillow_format)
                print(
                    "⚠️  SKIP existing: "
                    + path_image.name
                    + " | valid "
                    + pillow_format
                    + " "
                    + str(image_size[0])
                    + "x"
                    + str(image_size[1])
                )
                lpath_image.append(path_image)
                continue

            slide = presentation.Slides.Item(slide_number)
            slide.Export(str(path_image), powerpoint_format, width_px, height_px)
            original_size, cropped_size = _crop_image_whitespace(
                path_image,
                dpi,
                padding_mm,
                crop_threshold,
                pillow_format,
            )
            checked_size = _check_exported_image(path_image, pillow_format)
            if checked_size != cropped_size:
                raise RuntimeError(
                    "Image size changed during validation: " + path_image.name
                )
            print(
                "✅ "
                + path_ppt.name
                + " | slide "
                + str(slide_number)
                + " | "
                + str(original_size[0])
                + "x"
                + str(original_size[1])
                + " → "
                + str(cropped_size[0])
                + "x"
                + str(cropped_size[1])
            )
            lpath_image.append(path_image)
    except Exception as error:
        raise RuntimeError(
            "Failed to export " + str(path_ppt) + ": " + str(error)
        ) from error
    finally:
        if presentation is not None:
            presentation.Close()
        if powerpoint is not None:
            powerpoint.Quit()
        gc.collect()
        pythoncom.CoUninitialize()

    return lpath_image
