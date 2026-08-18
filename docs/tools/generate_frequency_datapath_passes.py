"""Generate self-contained weighted views of the frequency-calibration datapath.

The editable SVG is the single geometry source.  Visible SVG elements carry a
``data-semantic`` attribute; each reading pass assigns opacity weights to those
semantic groups.  Moving or restyling elements in Inkscape therefore propagates
to every generated pass without maintaining coordinate-based clipping masks.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


SVG_NS = "http://www.w3.org/2000/svg"
INKSCAPE_NS = "http://www.inkscape.org/namespaces/inkscape"
SODIPODI_NS = "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"

ET.register_namespace("", SVG_NS)
ET.register_namespace("inkscape", INKSCAPE_NS)
ET.register_namespace("sodipodi", SODIPODI_NS)

VISIBLE_TAGS = {"circle", "ellipse", "line", "path", "polygon", "polyline", "rect", "text"}
DEFAULT_WEIGHT = 0.14
CONTEXT_WEIGHT = 0.32


CURRENT_ID_SEMANTICS = {
    # Title, legend, lane headings, and footer.
    **{element_id: "chrome" for element_id in (
        "text4", "text5", "line5", "text6", "line6", "text7", "line7", "text8",
        "line8", "text9", "text10", "text19", "text47", "text84", "line124", "text125",
    )},
    # Control plane.
    **{element_id: "control-unit" for element_id in (
        "rect10", "text11", "text12", "text13", "rect13", "text14",
    )},
    **{element_id: "exception-unit" for element_id in (
        "rect14", "text15", "text16", "text17", "text18",
    )},
    # Register file and its independently weighted cells.
    **{element_id: "register-frame" for element_id in ("rect19", "rect20", "text20")},
    **{element_id: "register-sr" for element_id in ("rect21", "text21", "text22", "text23", "text24")},
    **{element_id: "register-er" for element_id in ("rect24", "text25", "text26", "text27", "text28")},
    **{element_id: "register-tr" for element_id in ("rect28", "text29", "text30", "text31", "text32")},
    **{element_id: "register-br" for element_id in ("rect32", "text33", "text34", "text35", "text36")},
    **{element_id: "register-cr" for element_id in ("rect36", "text37", "text38", "text39", "text40")},
    **{element_id: "register-ir" for element_id in ("rect40", "text41", "text42", "text43", "text44")},
    **{element_id: "snapshot-boundary" for element_id in ("rect44", "text45", "text46")},
    # Main command-event datapath.
    **{element_id: "operand-mux" for element_id in ("path47", "text48", "text49", "text50", "text51")},
    **{element_id: "cmd-register" for element_id in ("rect51", "text52", "text53", "text54")},
    **{element_id: "issue-unit" for element_id in ("rect54", "text55", "text56", "text57", "text58", "text59")},
    **{element_id: "measurement-unit" for element_id in ("rect59", "text60", "text61", "text62", "text63", "text64")},
    **{element_id: "evt-register" for element_id in ("rect64", "text65", "text66", "text67")},
    **{element_id: "identity-unit" for element_id in ("rect67", "text68", "text69", "text70", "text71", "text72")},
    **{element_id: "guard-comparator" for element_id in ("rect72", "text73", "text74", "text75", "text76", "text77")},
    **{element_id: "route-operands" for element_id in ("path77",)},
    **{element_id: "route-command-build" for element_id in ("path78",)},
    **{element_id: "route-command" for element_id in ("path79", "text79", "path80")},
    **{element_id: "route-event" for element_id in ("path81", "text81", "path82")},
    **{element_id: "route-id-valid" for element_id in ("path83", "text83")},
    # Update and persistence units.
    **{element_id: "secant-alu" for element_id in ("path84", "text85", "line85", "text86", "text87", "text88", "text89")},
    **{element_id: "commit-unit" for element_id in ("rect89", "text90", "text91", "text92", "text93", "text94", "text95")},
    **{element_id: "journal-memory" for element_id in ("rect95", "text96", "text97", "text98", "text99", "text100", "text101")},
    **{element_id: "checkpoint-memory" for element_id in ("rect101", "text102", "text103", "text104", "text105", "text106", "text107")},
    # Cross-unit signal routes.
    **{element_id: "route-context" for element_id in ("path107", "text108")},
    **{element_id: "route-lifecycle" for element_id in ("path108", "text109")},
    **{element_id: "route-control-bus" for element_id in ("path109", "path110", "text110")},
    **{element_id: "route-control-mux" for element_id in ("path111",)},
    **{element_id: "route-control-issue" for element_id in ("path112",)},
    **{element_id: "route-control-identity" for element_id in ("path113",)},
    **{element_id: "route-control-guard" for element_id in ("path114",)},
    **{element_id: "route-condition-codes" for element_id in ("path115", "text115")},
    **{element_id: "route-tracker-operands" for element_id in ("path116",)},
    **{element_id: "route-proposal" for element_id in ("path117", "text117")},
    **{element_id: "route-classified-event" for element_id in ("path118", "text118")},
    **{element_id: "route-tracker-event" for element_id in ("path119",)},
    **{element_id: "route-tr-writeback" for element_id in ("path120", "text120")},
    **{element_id: "route-state-writeback" for element_id in ("path121", "text121")},
    **{element_id: "route-journal" for element_id in ("path122", "text122")},
    **{element_id: "route-persist" for element_id in ("path123", "text123")},
    **{element_id: "route-restore" for element_id in ("path124", "text124")},
}


PASS_SPECS = {
    1: {
        "title": "PASS 1 · LANGUAGE & INTERFACES",
        "color": "#2b6ca3",
        "active": {
            "register-frame", "register-sr", "register-er", "register-tr", "register-br",
            "register-cr", "register-ir", "snapshot-boundary", "control-unit", "cmd-register",
            "evt-register", "identity-unit", "commit-unit", "route-context", "route-command-build",
            "route-command", "route-event", "route-id-valid", "route-state-writeback",
        },
        "context": {"operand-mux", "issue-unit", "measurement-unit", "guard-comparator"},
    },
    2: {
        "title": "PASS 2 · NORMAL PROTOCOL LOOP",
        "color": "#2b6ca3",
        "active": {
            "register-frame", "register-sr", "register-er", "register-br", "register-ir",
            "control-unit", "operand-mux", "cmd-register", "issue-unit", "measurement-unit",
            "evt-register", "identity-unit", "guard-comparator", "commit-unit", "route-context",
            "route-control-bus", "route-control-mux", "route-control-issue", "route-control-identity",
            "route-control-guard", "route-condition-codes", "route-operands", "route-command-build",
            "route-command", "route-event", "route-id-valid", "route-classified-event",
            "route-state-writeback",
        },
        "context": {"register-cr", "secant-alu", "journal-memory"},
    },
    3: {
        "title": "PASS 3 · MEASUREMENT BACKEND",
        "color": "#2b6ca3",
        "active": {
            "cmd-register", "issue-unit", "measurement-unit", "evt-register", "identity-unit",
            "route-command", "route-event", "route-id-valid", "route-control-issue",
            "route-control-identity",
        },
        "context": {"operand-mux", "guard-comparator", "control-unit", "commit-unit"},
    },
    4: {
        "title": "PASS 4 · TRACK FEEDBACK",
        "color": "#8a6d32",
        "active": {
            "register-tr", "register-cr", "register-ir", "operand-mux", "cmd-register", "issue-unit",
            "measurement-unit", "evt-register", "identity-unit", "guard-comparator", "secant-alu",
            "route-operands", "route-command-build", "route-command", "route-event", "route-id-valid",
            "route-control-mux", "route-control-issue", "route-control-identity", "route-control-guard",
            "route-tracker-operands", "route-proposal", "route-classified-event", "route-tracker-event",
            "route-tr-writeback",
        },
        "context": {"register-frame", "register-er", "control-unit", "commit-unit", "route-state-writeback"},
    },
    5: {
        "title": "PASS 5 · RUNTIME & PERSISTENCE",
        "color": "#40846f",
        "active": {
            "register-frame", "register-br", "register-ir", "snapshot-boundary", "control-unit",
            "issue-unit", "commit-unit", "journal-memory", "checkpoint-memory", "route-context",
            "route-control-issue", "route-state-writeback", "route-journal", "route-persist", "route-restore",
        },
        "context": {"cmd-register", "evt-register", "identity-unit", "guard-comparator"},
    },
    6: {
        "title": "PASS 6 · RECOVERY & EVIDENCE",
        "color": "#b0443c",
        "active": {
            "register-ir", "control-unit", "exception-unit", "operand-mux", "cmd-register", "issue-unit",
            "measurement-unit", "evt-register", "identity-unit", "guard-comparator", "commit-unit",
            "journal-memory", "checkpoint-memory", "route-lifecycle", "route-control-bus",
            "route-control-mux", "route-control-issue", "route-control-identity", "route-control-guard",
            "route-condition-codes", "route-operands", "route-command-build", "route-command", "route-event",
            "route-id-valid", "route-classified-event", "route-state-writeback", "route-journal",
            "route-persist", "route-restore",
        },
        "context": {"register-frame", "register-sr", "register-br", "snapshot-boundary"},
    },
}


def annotate_source(path: Path) -> None:
    source = path.read_text(encoding="utf-8")
    missing_ids = []
    for element_id, semantic in CURRENT_ID_SEMANTICS.items():
        tag_pattern = re.compile(
            rf'(<(?:text|rect|line|path|circle|ellipse|polygon|polyline)\b[^>]*\bid="{re.escape(element_id)}")([^>]*>)',
            re.DOTALL,
        )
        match = tag_pattern.search(source)
        if match is None:
            missing_ids.append(element_id)
            continue
        opening = match.group(0)
        if "data-semantic=" in opening:
            replacement = re.sub(
                r'data-semantic="[^"]*"',
                f'data-semantic="{semantic}"',
                opening,
            )
        else:
            replacement = f'{match.group(1)}\n       data-semantic="{semantic}"{match.group(2)}'
        source = source[: match.start()] + replacement + source[match.end() :]
    if missing_ids:
        raise ValueError("Cannot annotate missing SVG IDs: " + ", ".join(missing_ids))
    path.write_text(source, encoding="utf-8")


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _visible_elements(root: ET.Element) -> list[ET.Element]:
    parent_by_child = {child: parent for parent in root.iter() for child in parent}
    visible = []
    for element in root.iter():
        if _local_name(element.tag) not in VISIBLE_TAGS:
            continue
        parent = parent_by_child.get(element)
        inside_defs = False
        while parent is not None:
            if _local_name(parent.tag) in {"defs", "namedview"}:
                inside_defs = True
                break
            parent = parent_by_child.get(parent)
        if not inside_defs:
            visible.append(element)
    return visible


def validate_semantics(root: ET.Element) -> set[str]:
    visible = _visible_elements(root)
    untagged = [element.get("id", f"<{_local_name(element.tag)}>") for element in visible if not element.get("data-semantic")]
    if untagged:
        raise ValueError("Visible SVG elements lack data-semantic: " + ", ".join(untagged))
    semantics = {element.get("data-semantic") for element in visible}
    referenced = set().union(*(spec["active"] | spec["context"] for spec in PASS_SPECS.values()))
    missing = referenced - semantics
    if missing:
        raise ValueError("Pass configuration references missing semantics: " + ", ".join(sorted(missing)))
    return semantics


def _weight_for(semantic: str, spec: dict) -> float:
    if semantic in spec["active"]:
        return 1.0
    if semantic in spec["context"]:
        return CONTEXT_WEIGHT
    return DEFAULT_WEIGHT


def _append_pass_label(root: ET.Element, number: int, title: str, color: str) -> None:
    group = ET.SubElement(root, f"{{{SVG_NS}}}g", {
        "id": f"pass-{number}-label",
        f"{{{INKSCAPE_NS}}}groupmode": "layer",
        f"{{{INKSCAPE_NS}}}label": f"Pass {number} label",
    })
    ET.SubElement(group, f"{{{SVG_NS}}}rect", {
        "x": "1570", "y": "82", "width": "270", "height": "22", "rx": "11", "fill": color,
    })
    text = ET.SubElement(group, f"{{{SVG_NS}}}text", {
        "x": "1705", "y": "97", "text-anchor": "middle",
        "font-family": "Arial,Microsoft YaHei,sans-serif", "font-size": "12",
        "font-weight": "700", "style": "fill:#ffffff",
    })
    text.text = title


def _apply_weights(root: ET.Element, spec: dict) -> None:
    visible = _visible_elements(root)
    parent_by_child = {child: parent for parent in root.iter() for child in parent}
    for element in visible:
        parent = parent_by_child[element]
        semantic = element.get("data-semantic")
        wrapper = ET.Element(f"{{{SVG_NS}}}g", {
            "opacity": f"{_weight_for(semantic, spec):.2f}",
            "data-view-semantic": semantic,
        })
        element_index = list(parent).index(element)
        parent.remove(element)
        wrapper.append(element)
        parent.insert(element_index, wrapper)


def generate(source: Path, output_dir: Path) -> list[Path]:
    tree = ET.parse(source)
    root = tree.getroot()
    validate_semantics(root)
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = []
    for number, spec in PASS_SPECS.items():
        pass_tree = ET.ElementTree(ET.fromstring(ET.tostring(root, encoding="unicode")))
        pass_root = pass_tree.getroot()
        _apply_weights(pass_root, spec)
        _append_pass_label(pass_root, number, spec["title"], spec["color"])
        destination = output_dir / f"frequency_state_machine_pass{number}.svg"
        serialised = ET.tostring(pass_root, encoding="unicode")
        clean_lines = [line.rstrip() for line in serialised.splitlines()]
        destination.write_text(
            '<?xml version="1.0" encoding="utf-8"?>\n'
            + "\n".join(clean_lines)
            + "\n",
            encoding="utf-8",
        )
        generated.append(destination)
    return generated


def export_png(svg_paths: list[Path], inkscape: Path, width: int) -> None:
    for svg_path in svg_paths:
        destination = svg_path.with_suffix(".png")
        subprocess.run(
            [
                str(inkscape), str(svg_path), "--export-type=png",
                f"--export-filename={destination}", f"--export-width={width}",
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--annotate-source", action="store_true")
    parser.add_argument("--inkscape", type=Path)
    parser.add_argument("--png-width", type=int, default=1900)
    args = parser.parse_args()

    if args.annotate_source:
        annotate_source(args.source)
    generated = generate(args.source, args.output_dir)
    if args.inkscape:
        export_png(generated, args.inkscape, args.png_width)


if __name__ == "__main__":
    main()
