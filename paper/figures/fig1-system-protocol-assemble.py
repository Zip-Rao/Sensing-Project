#!/usr/bin/env python3
"""Assemble the four editable Fig. 1 SVG panels into one publication SVG."""

from __future__ import annotations

import copy
import re
import xml.etree.ElementTree as ET
from pathlib import Path


FIGURES = Path(r"C:\Users\21034\Desktop\Workspace\Sensing project\Sensing-Project\paper\figures")
OUT = Path(r"C:\Users\21034\.codex\visualizations\2026\08\08\019fe034-5301-7a33-b8ba-b7d6795c75e5")

SVG = "http://www.w3.org/2000/svg"
INK = "http://www.inkscape.org/namespaces/inkscape"
XLINK = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG)
ET.register_namespace("inkscape", INK)
ET.register_namespace("xlink", XLINK)


PANELS = [
    # file, layer name, x, y, displayed width in mm
    ("fig1a-v2-editable.svg", "Panel a - global dispersion", 0.0, 0.0, 90.0),
    ("fig1b-v6-editable.svg", "Panel b - probe timing", 94.0, 0.0, 78.0),
    ("fig1c-v1-editable.svg", "Panel c - differential response", 0.0, 61.3, 90.0),
    ("fig1d-v8-editable.svg", "Panel d - calibration workflow", 94.0, 61.3, 86.0),
]

CANVAS_W = 180.0
CANVAS_H = 106.3


def parse_viewbox(root: ET.Element) -> tuple[float, float, float, float]:
    raw = root.attrib.get("viewBox")
    if not raw:
        raise ValueError("source SVG has no viewBox")
    values = tuple(float(v) for v in re.split(r"[ ,]+", raw.strip()))
    if len(values) != 4:
        raise ValueError(f"invalid viewBox: {raw}")
    return values  # type: ignore[return-value]


def prefix_ids(node: ET.Element, prefix: str) -> None:
    id_map: dict[str, str] = {}
    for elem in node.iter():
        old = elem.attrib.get("id")
        if old:
            new = f"{prefix}-{old}"
            id_map[old] = new
            elem.set("id", new)

    def replace_refs(value: str) -> str:
        for old, new in id_map.items():
            value = value.replace(f"url(#{old})", f"url(#{new})")
            if value == f"#{old}":
                value = f"#{new}"
        return value

    for elem in node.iter():
        for key, value in list(elem.attrib.items()):
            elem.set(key, replace_refs(value))
        if elem.text and "url(#" in elem.text:
            elem.text = replace_refs(elem.text)


def remove_source_panel_labels(node: ET.Element) -> None:
    labels = {"(a)", "(b)", "(c)", "(d)"}
    for parent in node.iter():
        for child in list(parent):
            text = "".join(child.itertext()).strip()
            if child.tag == f"{{{SVG}}}text" and text in labels:
                parent.remove(child)


def add_panel(composite: ET.Element, filename: str, label: str, x: float, y: float, width: float, index: int) -> float:
    source_root = ET.parse(FIGURES / filename).getroot()
    min_x, min_y, vb_w, vb_h = parse_viewbox(source_root)
    scale = width / vb_w
    height = vb_h * scale

    layer = ET.SubElement(
        composite,
        f"{{{SVG}}}g",
        {
            "id": f"panel-{index}",
            f"{{{INK}}}groupmode": "layer",
            f"{{{INK}}}label": label,
            "transform": f"translate({x:.6f},{y:.6f}) scale({scale:.9f}) translate({-min_x:.6f},{-min_y:.6f})",
        },
    )
    for child in list(source_root):
        if child.tag in {f"{{{SVG}}}title", f"{{{SVG}}}desc", f"{{{SVG}}}metadata"}:
            continue
        layer.append(copy.deepcopy(child))
    remove_source_panel_labels(layer)
    prefix_ids(layer, f"p{index}")
    return height


def main() -> None:
    root = ET.Element(
        f"{{{SVG}}}svg",
        {
            "width": f"{CANVAS_W}mm",
            "height": f"{CANVAS_H}mm",
            "viewBox": f"0 0 {CANVAS_W} {CANVAS_H}",
            "version": "1.1",
        },
    )
    ET.SubElement(root, f"{{{SVG}}}title").text = "Fig. 1: system, probes, response, and calibration workflow"
    ET.SubElement(root, f"{{{SVG}}}desc").text = (
        "Editable two-column assembly of panels (a)-(d). Panel groups are separate Inkscape layers."
    )
    background = ET.SubElement(
        root,
        f"{{{SVG}}}g",
        {
            "id": "background",
            f"{{{INK}}}groupmode": "layer",
            f"{{{INK}}}label": "Background",
        },
    )
    ET.SubElement(
        background,
        f"{{{SVG}}}rect",
        {"x": "0", "y": "0", "width": str(CANVAS_W), "height": str(CANVAS_H), "fill": "#ffffff"},
    )

    heights = []
    for index, (filename, label, x, y, width) in enumerate(PANELS, start=1):
        heights.append((filename, add_panel(root, filename, label, x, y, width, index)))

    labels = ET.SubElement(
        root,
        f"{{{SVG}}}g",
        {
            "id": "panel-labels",
            f"{{{INK}}}groupmode": "layer",
            f"{{{INK}}}label": "Panel labels",
            "font-family": '"Times New Roman", "Latin Modern Roman", serif',
            "font-size": "3.2",
            "font-weight": "700",
            "fill": "#17191c",
        },
    )
    for text, x, y in (("(a)", 0.8, 3.4), ("(b)", 94.8, 3.4), ("(c)", 0.8, 64.7), ("(d)", 94.8, 64.7)):
        ET.SubElement(labels, f"{{{SVG}}}text", {"x": str(x), "y": str(y)}).text = text

    out = OUT / "fig1-system-protocol-v1-editable.svg"
    ET.ElementTree(root).write(out, encoding="utf-8", xml_declaration=True)
    print(f"wrote {out}")
    for filename, height in heights:
        print(f"  {filename}: displayed height {height:.2f} mm")


if __name__ == "__main__":
    main()
