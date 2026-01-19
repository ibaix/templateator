"""JSON template reader for ring mesh generation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Tuple

from .ring_template_reader import (
    NodeRecord,
    ElementRecord,
    MeshTemplate,
)


class JsonTemplateError(RuntimeError):
    """Exception for JSON template format errors."""


# Element type definitions: type -> (expected_nodes, order, 3D_type)
ELEMENT_TYPES = {
    # First-order (linear) elements
    "QUAD": (4, 1, "hex8"),
    "QUAD4": (4, 1, "hex8"),
    "TRI": (3, 1, "penta6"),
    "TRI3": (3, 1, "penta6"),
    # Second-order (quadratic) elements
    "QUAD8": (8, 2, "hex20"),
    "TRI6": (6, 2, "penta15"),
}


def read_json(path: str | Path) -> MeshTemplate:
    """Read a JSON ring template file.

    Parameters
    ----------
    path:
        Path to the JSON template file.

    Returns
    -------
    MeshTemplate
        Parsed template with nodes, elements, and surfaces.
    """

    filepath = Path(path)
    if not filepath.exists():
        raise FileNotFoundError(f"Template not found: {filepath}")

    with filepath.open("r", encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            raise JsonTemplateError(f"Invalid JSON: {e}") from e

    # Validate required fields
    required_fields = ["name", "version", "nodes", "elements"]
    for field in required_fields:
        if field not in data:
            raise JsonTemplateError(f"Missing required field: {field}")

    # Parse nodes
    nodes: Dict[int, NodeRecord] = {}
    for node_data in data["nodes"]:
        node_id = int(node_data["id"])
        if node_id in nodes:
            raise JsonTemplateError(f"Duplicate node ID: {node_id}")
        nodes[node_id] = NodeRecord(
            id=node_id,
            x=float(node_data["x"]),
            y=float(node_data.get("y", 0.0)),
            z=float(node_data["z"]),
        )

    # Parse elements
    elements: Dict[int, ElementRecord] = {}
    element_type: str | None = None
    element_order: int = 1  # Track element order (1=linear, 2=quadratic)

    for elem_data in data["elements"]:
        elem_id = int(elem_data["id"])
        if elem_id in elements:
            raise JsonTemplateError(f"Duplicate element ID: {elem_id}")

        elem_type = elem_data["type"].upper()
        connectivity = tuple(int(n) for n in elem_data["nodes"])

        # Look up element type definition
        if elem_type not in ELEMENT_TYPES:
            raise JsonTemplateError(
                f"Element {elem_id}: Unknown element type '{elem_type}'. "
                f"Supported types: {', '.join(ELEMENT_TYPES.keys())}"
            )

        expected_nodes, order, febio_type = ELEMENT_TYPES[elem_type]

        # Validate connectivity length
        if len(connectivity) != expected_nodes:
            raise JsonTemplateError(
                f"Element {elem_id}: {elem_type} requires {expected_nodes} nodes, "
                f"got {len(connectivity)}"
            )

        # Validate node references
        for nid in connectivity:
            if nid not in nodes:
                raise JsonTemplateError(
                    f"Element {elem_id} references unknown node {nid}"
                )

        # Set 3D element type (what the shell becomes after revolution)
        if element_type is None:
            element_type = febio_type
            element_order = order

        elements[elem_id] = ElementRecord(
            id=elem_id,
            connectivity=connectivity,
            elset=data.get("name", "ring_2D"),
        )

    if not nodes:
        raise JsonTemplateError("Template contains no nodes")
    if not elements:
        raise JsonTemplateError("Template contains no elements")

    # Parse surfaces (optional)
    surfaces: Dict[str, Tuple[Tuple[int, int], ...]] = {}
    if "surfaces" in data:
        for surface_name, entries in data["surfaces"].items():
            surface_entries = []
            for entry in entries:
                elem_id = int(entry["element"])
                face_idx = int(entry["face"])
                surface_entries.append((elem_id, face_idx))
            surfaces[surface_name] = tuple(surface_entries)

    return MeshTemplate(
        nodes=nodes,
        elements=elements,
        element_type=element_type or "hex8",
        surfaces=surfaces,
        element_order=element_order,
    )
