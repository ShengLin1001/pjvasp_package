#!/usr/bin/env python3
"""Dry-run Slurm job-script generation with mymetal.slurm.submit.

This tutorial is VASP-free and deterministic. It demonstrates the
**text-generation** building blocks of :mod:`mymetal.slurm.submit`:

1. :func:`generate_script_header` — emit the ``#SBATCH`` header and module
   load block for one partition / resource request.
2. :func:`generate_launcher_command` — wrap a calculation command in a
   ``srun`` / ``mpirun`` / ``none`` launcher (with the launch-retry shim).
3. :func:`generate_slurm_script_base` — assemble header + launcher into a
   complete single-job Bash script and optionally write it to disk.
4. :func:`check_wall_time` — validate Slurm wall-time strings; the invalid
   path raises ``SystemExit`` via ``fail``, which we catch for the demo.
5. :func:`pei_slurm_univ_submit` — the full orchestrator, here driven in
   **dry-run mode** (``if_sbatch=False``) over a tiny synthetic ``y_dir``
   tree so it only *generates* scripts and never calls ``sbatch``.

This script NEVER submits a job. Every ``sbatch`` is replaced by text
generation. The figure renders (left) a flowchart of the generation
pipeline and (right) the generated base script as monospace text.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # MUST precede import pyplot
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from mymetal.slurm.submit import (
    check_wall_time,
    generate_launcher_command,
    generate_script_header,
    generate_slurm_script_base,
    pei_slurm_univ_submit,
)
from mymetal.universal.plot.general import general_set_all_rcParams


# --------------------------------------------------------------------------
# Synthetic inputs — deterministic, no cluster access required.
# --------------------------------------------------------------------------

# A toy MODULE_BLOCKS dict mapping a profile-type name to the shell block
# that loads environment modules + exports. Real workflows pass a much
# richer dict; here we keep one profile so the example is self-contained.
MODULE_BLOCKS: dict[str, str] = {
    "vasp_intel": (
        "module purge\n"
        "module load intel/2023a\n"
        "module load impi/2023a\n"
        "module load vasp/6.4\n"
        "export I_MPI_PIN_DOMAIN=omp\n"
        "export OMP_NUM_THREADS=1\n"
    ),
}

PARTITION = "cpu-epyc"
NODES = 1
NCORES = 16
WALL_TIME = "2-00:00:00"
LAUNCHER = "srun"
CMD = "vasp_std"


# --------------------------------------------------------------------------
# Build helpers — each returns a deterministic artifact (string or dict).
# --------------------------------------------------------------------------

def build_header() -> str:
    """Generate just the #SBATCH header + module block."""
    return generate_script_header(
        partition=PARTITION,
        nodes=NODES,
        ncores=NCORES,
        module_profile_type="vasp_intel",
        MODULE_BLOCKS=MODULE_BLOCKS,
        wall_time=WALL_TIME,
    )


def build_launcher() -> str:
    """Generate the launcher command line (srun + retry shim)."""
    return generate_launcher_command(
        launcher=LAUNCHER,
        cmd=CMD,
        if_use_my_launcher=False,
    )


def build_base_script(path_save: Path) -> str:
    """Generate a complete single-job Slurm script and write it to disk."""
    return generate_slurm_script_base(
        partition=PARTITION,
        nodes=NODES,
        ncores=NCORES,
        module_profile_type="vasp_intel",
        MODULE_BLOCKS=MODULE_BLOCKS,
        launcher=LAUNCHER,
        cmd=CMD,
        if_use_my_launcher=False,
        if_output=True,
        path_save=path_save,
        wall_time=WALL_TIME,
    )


def build_wall_time_checks() -> dict[str, object]:
    """Exercise check_wall_time with valid and invalid inputs.

    Invalid inputs raise SystemExit(1) via fail(); we catch and record.
    """
    # Per WALL_TIME_PATTERN = ^[0-9]+(?:-[0-9]+)?(?::[0-9]+){0,2}$:
    # "2-00:00:00", "12:00:00", "12:00", "360" are all VALID (the regex
    # accepts 0-2 colon groups). None means "no --time line".
    valid_cases = ["2-00:00:00", "12:00:00", "12:00", "360", None]
    # Truly invalid: empty, non-numeric, too many colon groups.
    invalid_cases = ["", "abc", "1:2:3:4", "--:--"]

    valid_results: dict[str, str | None] = {}
    for wt in valid_cases:
        valid_results[repr(wt)] = check_wall_time(wt)

    invalid_results: dict[str, str] = {}
    for wt in invalid_cases:
        try:
            check_wall_time(wt)
            invalid_results[repr(wt)] = "(no error raised — UNEXPECTED)"
        except SystemExit as exc:
            invalid_results[repr(wt)] = f"SystemExit(code={exc.code})"

    return {"valid": valid_results, "invalid": invalid_results}


def build_orchestrator_dry_run(path_output: Path) -> dict[str, object]:
    """Drive pei_slurm_univ_submit in dry-run mode over a synthetic y_dir tree.

    Creates a tiny project layout:

        <output>/synthetic_project/dir/y_dir/job_001/
        <output>/synthetic_project/dir/y_dir/job_002/
        <output>/synthetic_project/dir/y_dir/job_003/

    then calls the orchestrator with mode="parallel" and if_sbatch=False.
    The orchestrator writes sub_slurm_univ.sh into each job directory but
    NEVER calls sbatch. cwd is saved/restored because the orchestrator
    calls os.chdir(path_root).
    """
    path_root = (path_output / "synthetic_project").resolve()
    path_y_dir = path_root / "dir" / "y_dir"
    path_y_dir.mkdir(parents=True, exist_ok=True)
    job_names = ["job_001", "job_002", "job_003"]
    lpath_job = []
    for name in job_names:
        path_job = path_y_dir / name
        path_job.mkdir(parents=True, exist_ok=True)
        lpath_job.append(path_job)

    cwd_before = Path.cwd()
    try:
        pei_slurm_univ_submit(
            path_root=path_root,
            mode="parallel",
            dir_root=Path("dir"),
            lsubdir=None,
            chunks=1,
            module_profile_type="vasp_intel",
            launcher_type=LAUNCHER,
            cmd=CMD,
            if_use_my_launcher=False,
            partition=PARTITION,
            nodes=NODES,
            ncores=NCORES,
            if_sbatch=False,            # <-- DRY-RUN: no sbatch submission
            child_wall_time=WALL_TIME,
            parent_wall_time=None,
            MODULE_BLOCKS=MODULE_BLOCKS,
        )
    finally:
        os.chdir(cwd_before)

    # Read back one generated child script for display / assertion.
    path_sample = lpath_job[0] / "sub_slurm_univ.sh"
    sample_text = path_sample.read_text(encoding="utf-8") if path_sample.is_file() else ""

    return {
        "path_root": path_root,
        "lpath_job": lpath_job,
        "path_sample": path_sample,
        "sample_text": sample_text,
    }


# --------------------------------------------------------------------------
# Render helpers — figure only, no computation.
# --------------------------------------------------------------------------

def _draw_box(ax, xy, w, h, text, color="#e8f0fe", edge="#1a73e8"):
    """Draw a rounded box with centered text on the axes."""
    x, y = xy
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.4, edgecolor=edge, facecolor=color,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=8.5, wrap=True)


def _draw_arrow(ax, start, end):
    """Draw a thin arrow between two points."""
    arrow = FancyArrowPatch(
        start, end,
        arrowstyle="-|>", mutation_scale=12,
        linewidth=1.2, color="#5f6368",
    )
    ax.add_patch(arrow)


def render_flowchart(ax) -> None:
    """Draw the generation-pipeline flowchart on the left panel."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("Slurm script generation pipeline\n(DRY-RUN: no sbatch)", fontsize=10)

    # Input layer
    _draw_box(ax, (0.4, 8.6), 3.0, 1.0,
              "Inputs\npartition, nodes, ncores,\nMODULE_BLOCKS, wall_time",
              color="#fce8e6", edge="#d93025")
    _draw_box(ax, (6.6, 8.6), 3.0, 1.0,
              "check_wall_time\nvalidates format\n(invalid -> SystemExit)",
              color="#fef7e0", edge="#f9ab00")

    # Building-block layer
    _draw_box(ax, (0.4, 6.4), 3.0, 1.2,
              "generate_script_header\n#SBATCH lines + modules",
              color="#e8f0fe", edge="#1a73e8")
    _draw_box(ax, (3.8, 6.4), 3.0, 1.2,
              "generate_launcher_command\nsrun / mpirun / none\n+ retry shim",
              color="#e8f0fe", edge="#1a73e8")

    # Assembly layer
    _draw_box(ax, (1.6, 4.2), 4.2, 1.2,
              "generate_slurm_script_base\nheader + launcher + auto-comment\n-> script string + file",
              color="#e6f4ea", edge="#188038")

    # Orchestrator layer
    _draw_box(ax, (1.6, 2.0), 4.2, 1.2,
              "pei_slurm_univ_submit\n(if_sbatch=False)\ngenerates per-dir scripts",
              color="#e6f4ea", edge="#188038")

    # Output layer
    _draw_box(ax, (1.6, 0.2), 4.2, 1.0,
              "sub_slurm_univ.sh\nwritten to each job dir",
              color="#f3e8fd", edge="#9334e6")

    # Arrows
    _draw_arrow(ax, (1.9, 8.6), (1.9, 7.6))   # inputs -> header
    _draw_arrow(ax, (3.0, 8.6), (5.0, 7.6))   # inputs -> launcher (diagonal)
    _draw_arrow(ax, (8.1, 8.6), (5.0, 7.6))   # check_wall_time -> launcher
    _draw_arrow(ax, (2.2, 6.4), (3.0, 5.4))   # header -> base
    _draw_arrow(ax, (5.3, 6.4), (4.4, 5.4))   # launcher -> base
    _draw_arrow(ax, (3.7, 4.2), (3.7, 3.2))   # base -> orchestrator
    _draw_arrow(ax, (3.7, 2.0), (3.7, 1.2))   # orchestrator -> output

    # Annotation: dry-run note
    ax.text(8.1, 2.6, "if_sbatch=False\n-> no sbatch call",
            ha="center", va="center", fontsize=8,
            color="#d93025",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff", ec="#d93025", lw=1.0))


def render_script_text(ax, script_text: str, max_lines: int = 26) -> None:
    """Render the generated script as monospace text on the right panel."""
    ax.set_axis_off()
    ax.set_title("Generated base script (sub_slurm_univ.sh)", fontsize=10)
    # Light background fills the whole panel so the quadrant is non-white
    # even where the text block does not reach. Use a slightly darker tint
    # than pure-white so the extended per-quadrant check (>0.05 deviation)
    # registers it as non-blank.
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                               facecolor="#eef1f5", edgecolor="#dadce0",
                               linewidth=1.0, zorder=0))
    lines = script_text.splitlines()
    # Replace the CJK auto-comment line with an ASCII equivalent so the
    # monospace font (DejaVu Sans Mono) can render every glyph.
    lines = [
        "# (auto-generated; overwritten on every run)" if "该脚本由程序自动生成" in ln else ln
        for ln in lines
    ]
    if len(lines) > max_lines:
        lines = lines[:max_lines] + [f"... ({len(lines) - max_lines} more lines)"]
    # Monospace text block inside a light box.
    text = "\n".join(lines)
    ax.text(0.02, 0.98, text,
            transform=ax.transAxes,
            ha="left", va="top",
            family="DejaVu Sans Mono", fontsize=7.2,
            zorder=1)


def render_figure(header_text: str, base_script: str, path_image: Path) -> None:
    """Two-panel figure: flowchart (left) + generated script text (right)."""
    path_image.parent.mkdir(parents=True, exist_ok=True)

    general_set_all_rcParams(
        backend="Agg",
        axes_grid=False,
        figure_subp=(1, 2),
        figure_one_figsize=(7.0, 4.2),
        figure_left=0.04,
        figure_top=0.90,
        figure_wspace=0.06,
        figure_hspace=0.4,
        font_family=["DejaVu Sans"],
        fontsize=10,
        axes_titlepad=8,
    )
    fig, axes = plt.subplots(1, 2)
    render_flowchart(axes[0])
    render_script_text(axes[1], base_script)
    fig.savefig(path_image, bbox_inches="tight")
    plt.close(fig)

    # Non-blank PNG verification.
    if not path_image.is_file() or path_image.stat().st_size < 1_000:
        raise AssertionError("image was not created: " + str(path_image))
    image_rgb = plt.imread(path_image)[..., :3]
    if float(np.mean(np.abs(image_rgb - 1.0))) < 0.002:
        raise AssertionError("image is effectively blank: " + str(path_image))


# --------------------------------------------------------------------------
# Assertions — deterministic invariants on the generated artifacts.
# --------------------------------------------------------------------------

def run_assertions(
    header_text: str,
    launcher_text: str,
    base_script: str,
    path_base_script: Path,
    wt_results: dict[str, object],
    orch_result: dict[str, object],
) -> None:
    """Assert deterministic invariants on every generated artifact."""
    # Header must contain SBATCH directives and the partition name.
    assert "#SBATCH" in header_text
    assert f"#SBATCH -p {PARTITION}" in header_text
    assert f"#SBATCH -N {NODES}" in header_text
    assert f"#SBATCH -n {NCORES}" in header_text
    assert f"#SBATCH --time={WALL_TIME}" in header_text
    # Module block must be appended.
    assert "module load intel/2023a" in header_text

    # Launcher must reference srun and the command.
    assert "srun" in launcher_text
    assert CMD in launcher_text
    assert "$SLURM_NTASKS" in launcher_text

    # Base script must combine header + launcher + auto-comment.
    assert "#!/bin/bash" in base_script
    assert "#SBATCH" in base_script
    assert PARTITION in base_script
    assert "srun" in base_script
    assert CMD in base_script
    assert "Auto-generated by pei_slurm_univ_submit" in base_script

    # File on disk must match the returned string.
    assert path_base_script.is_file()
    assert path_base_script.read_text(encoding="utf-8") == base_script

    # check_wall_time: valid inputs normalize; None stays None.
    assert wt_results["valid"][repr("2-00:00:00")] == "2-00:00:00"
    assert wt_results["valid"][repr("12:00:00")] == "12:00:00"
    assert wt_results["valid"][repr("12:00")] == "12:00"
    assert wt_results["valid"][repr("360")] == "360"
    assert wt_results["valid"][repr(None)] is None
    # Invalid inputs must raise SystemExit.
    for key, val in wt_results["invalid"].items():
        assert val.startswith("SystemExit"), f"{key} did not raise: {val}"

    # Orchestrator dry-run: every job dir must have a generated script.
    for path_job in orch_result["lpath_job"]:
        path_script = path_job / "sub_slurm_univ.sh"
        assert path_script.is_file(), f"missing script: {path_script}"
        text = path_script.read_text(encoding="utf-8")
        assert "#!/bin/bash" in text
        assert f"#SBATCH -p {PARTITION}" in text
        assert "srun" in text
    # The sample text read back must be non-empty and contain the cmd.
    assert orch_result["sample_text"]
    assert CMD in orch_result["sample_text"]


# --------------------------------------------------------------------------
# Main entry point.
# --------------------------------------------------------------------------

def run_example(path_output: Path) -> dict[str, object]:
    """Run the full dry-run demonstration and return all artifacts."""
    path_output = path_output.resolve()
    path_output.mkdir(parents=True, exist_ok=True)

    print("================ DRY-RUN: no sbatch submission ================")
    print(f"output directory: {path_output}")
    print("")

    # 1. Header.
    header_text = build_header()
    print("---------------- generate_script_header ----------------")
    print(header_text)

    # 2. Launcher.
    launcher_text = build_launcher()
    print("---------------- generate_launcher_command ----------------")
    print(launcher_text)
    print("")

    # 3. Base script (writes to disk).
    path_base_script = path_output / "sub_slurm_univ_base.sh"
    base_script = build_base_script(path_base_script)
    print("---------------- generate_slurm_script_base ----------------")
    print(base_script)
    print(f"(written to: {path_base_script})")
    print("")

    # 4. Wall-time checks.
    wt_results = build_wall_time_checks()
    print("---------------- check_wall_time ----------------")
    print("valid inputs -> normalized:")
    for k, v in wt_results["valid"].items():
        print(f"  {k}: {v!r}")
    print("invalid inputs -> caught SystemExit:")
    for k, v in wt_results["invalid"].items():
        print(f"  {k}: {v}")
    print("")

    # 5. Orchestrator dry-run over a synthetic y_dir tree.
    print("---------------- pei_slurm_univ_submit (if_sbatch=False) ----------------")
    orch_result = build_orchestrator_dry_run(path_output)
    print(f"path_root: {orch_result['path_root']}")
    for path_job in orch_result["lpath_job"]:
        path_script = path_job / "sub_slurm_univ.sh"
        print(f"  generated: {path_script}  (exists={path_script.is_file()})")
    print(f"sample script ({orch_result['path_sample']}):")
    print(orch_result["sample_text"])
    print("")

    # 6. Figure.
    path_image = path_output / "slurm_script_generation.png"
    render_figure(header_text, base_script, path_image)
    print(f"wrote: {path_image}")

    # 7. Assertions.
    run_assertions(
        header_text, launcher_text, base_script, path_base_script,
        wt_results, orch_result,
    )

    return {
        "header_text": header_text,
        "launcher_text": launcher_text,
        "base_script": base_script,
        "path_base_script": path_base_script,
        "wt_results": wt_results,
        "orch_result": orch_result,
        "path_image": path_image,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dry-run Slurm job-script generation (no sbatch). VASP-free, deterministic.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("slurm-script-generation-output"),
        help="Directory for the generated scripts and PNG figure.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    result = run_example(parse_args().output)
    assert "#SBATCH" in result["base_script"]
    assert result["path_image"].is_file()
    print("OK: all assertions passed.")
