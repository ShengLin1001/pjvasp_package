#!/usr/bin/env python3
"""Regenerate deterministic structure images used by the Sphinx manual."""

from __future__ import annotations

import sys
from pathlib import Path


PATH_DOCS = Path(__file__).resolve().parents[1]
PATH_EXAMPLES = PATH_DOCS / "examples"
PATH_IMAGES = PATH_DOCS / "source" / "_static" / "images" / "generated"
sys.path.insert(0, str(PATH_EXAMPLES))

from getting_started_au111 import build_au111_slab, render_slab


def main() -> None:
    PATH_IMAGES.mkdir(parents=True, exist_ok=True)
    path_image = PATH_IMAGES / "au111_slab.png"
    render_slab(build_au111_slab(), path_image)
    print("generated: " + str(path_image))


if __name__ == "__main__":
    main()
