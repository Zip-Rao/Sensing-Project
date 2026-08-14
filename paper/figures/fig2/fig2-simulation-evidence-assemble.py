#!/usr/bin/env python3
"""Assemble the four simulation-evidence panels into an editable 2x2 SVG."""
from __future__ import annotations

import copy
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


HERE = Path(__file__).resolve().parent
INKSCAPE = Path(r"C:\Program Files\Inkscape\bin\inkscape.exe")
SVG = "http://www.w3.org/2000/svg"
INK = "http://www.inkscape.org/namespaces/inkscape"
ET.register_namespace("", SVG)
ET.register_namespace("inkscape", INK)

PANELS = (
    ("fig2a-short-pulse-response.svg", "Panel a - response", 0.0, 0.0, 86.0),
    ("fig2b-frequency-error.svg", "Panel b - flux recovery", 94.0, 0.0, 86.0),
    ("fig2c-validation-summary.svg", "Panel c - validation", 0.0, 73.0, 86.0),
    ("fig2d-accuracy-cost.svg", "Panel d - solver cost", 94.0, 73.0, 86.0),
)
CANVAS_W = 180.0
CANVAS_H = 142.0


def viewbox(root: ET.Element) -> tuple[float, float, float, float]:
    values = tuple(float(value) for value in re.split(r"[ ,]+", root.attrib["viewBox"].strip()))
    if len(values) != 4:
        raise ValueError("invalid source SVG viewBox")
    return values  # type: ignore[return-value]


def prefix_ids(node: ET.Element, prefix: str) -> None:
    mapping = {}
    for element in node.iter():
        old = element.attrib.get("id")
        if old:
            mapping[old] = f"{prefix}-{old}"
            element.set("id", mapping[old])
    for element in node.iter():
        for key, value in list(element.attrib.items()):
            for old, new in mapping.items():
                value = value.replace(f"url(#{old})", f"url(#{new})")
                if value == f"#{old}":
                    value = f"#{new}"
            element.set(key, value)


def add_panel(root: ET.Element, filename: str, label: str, x: float, y: float,
              width: float, index: int) -> None:
    source = ET.parse(HERE / filename).getroot()
    min_x, min_y, vb_w, _ = viewbox(source)
    scale = width / vb_w
    layer = ET.SubElement(root, f"{{{SVG}}}g", {
        "id": f"panel-{index}",
        f"{{{INK}}}groupmode": "layer",
        f"{{{INK}}}label": label,
        "transform": f"translate({x:.6f},{y:.6f}) scale({scale:.9f}) "
                     f"translate({-min_x:.6f},{-min_y:.6f})",
    })
    for child in list(source):
        if child.tag not in {f"{{{SVG}}}title", f"{{{SVG}}}desc", f"{{{SVG}}}metadata"}:
            layer.append(copy.deepcopy(child))
    prefix_ids(layer, f"p{index}")


def export(svg_path: Path) -> None:
    if not INKSCAPE.exists():
        raise FileNotFoundError(f"Inkscape not found: {INKSCAPE}")
    for suffix, options in (
        ("pdf", []),
        ("png", ["--export-dpi=600"]),
    ):
        output = svg_path.with_suffix(f".{suffix}")
        subprocess.run([
            str(INKSCAPE), str(svg_path), f"--export-filename={output}", *options,
        ], check=True)


def main() -> None:
    root = ET.Element(f"{{{SVG}}}svg", {
        "width": f"{CANVAS_W}mm", "height": f"{CANVAS_H}mm",
        "viewBox": f"0 0 {CANVAS_W} {CANVAS_H}", "version": "1.1",
    })
    ET.SubElement(root, f"{{{SVG}}}title").text = "Simulation-only Fig. 2 evidence"
    ET.SubElement(root, f"{{{SVG}}}desc").text = (
        "Deterministic numerical evidence; not an experimental main-text figure."
    )
    background = ET.SubElement(root, f"{{{SVG}}}g", {
        "id": "background", f"{{{INK}}}groupmode": "layer",
        f"{{{INK}}}label": "Background",
    })
    ET.SubElement(background, f"{{{SVG}}}rect", {
        "x": "0", "y": "0", "width": str(CANVAS_W), "height": str(CANVAS_H),
        "fill": "#ffffff",
    })
    for index, panel in enumerate(PANELS, start=1):
        add_panel(root, *panel, index)
    output = HERE / "fig2-simulation-evidence.svg"
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    export(output)
    print(f"wrote {output.with_suffix('.*')}")


if __name__ == "__main__":
    main()
