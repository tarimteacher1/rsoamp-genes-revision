#!/usr/bin/env python3
"""Stack Figure 2A-C into one submission image while preserving vector PDF content."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf._page import PageObject


STEMS = ["Figure2A_Defensin_phylogeny", "Figure2B_Snakin_GASA_phylogeny", "Figure2C_nsLTP_phylogeny"]


def combine_png(folder: Path, output: Path):
    images = [Image.open(folder / f"{stem}.png").convert("RGB") for stem in STEMS]
    gap = 180
    width = max(image.width for image in images)
    height = sum(image.height for image in images) + gap * (len(images) - 1)
    canvas = Image.new("RGB", (width, height), "white")
    y = 0
    for image in images:
        canvas.paste(image, ((width - image.width) // 2, y))
        y += image.height + gap
    canvas.save(output, dpi=(600, 600), optimize=True)


def combine_pdf(folder: Path, output: Path):
    pages = [PdfReader(str(folder / f"{stem}.pdf")).pages[0] for stem in STEMS]
    gap = 18.0
    widths = [float(page.mediabox.width) for page in pages]
    heights = [float(page.mediabox.height) for page in pages]
    width = max(widths)
    height = sum(heights) + gap * (len(pages) - 1)
    combined = PageObject.create_blank_page(width=width, height=height)
    y_top = height
    for page, page_width, page_height in zip(pages, widths, heights):
        y_top -= page_height
        combined.merge_translated_page(page, (width - page_width) / 2.0, y_top, expand=False)
        y_top -= gap
    writer = PdfWriter()
    writer.add_page(combined)
    with output.open("wb") as handle:
        writer.write(handle)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    combine_png(args.figure_dir, args.figure_dir / "Figure2_combined_phylogenies.png")
    combine_pdf(args.figure_dir, args.figure_dir / "Figure2_combined_phylogenies.pdf")
    print("PASS Figure2 combined PNG/PDF")


if __name__ == "__main__":
    main()
