#Author-KIROSHI
#Description-Export sketch points as JSON ring template for KIROSHI mesh generator

# pyright: reportMissingImports=false
# type: ignore[import]

"""
Fusion 360 Add-In Script: Ring Template JSON Exporter

Run the script from within Fusion 360 (Scripts and Add-Ins > Shift+S).

This script extracts sketch points from a Fusion 360 sketch and exports them
as a JSON template file compatible with the KIROSHI mesh generator.

USAGE:
1. Open Fusion 360 with your ring cross-section sketch
2. Run this script from Scripts and Add-Ins (Shift+S)
3. Select the sketch containing the ring profile points
4. Define elements interactively or modify PREDEFINED_ELEMENTS below
5. Choose output file location

COORDINATE SYSTEM:
- Fusion 360 coordinate units are configurable (see FUSION_UNITS below)
- Fusion XY plane maps to template XZ plane (Y=0)
- X = radial direction (inner to outer)
- Z = thickness direction (bottom to top)

ELEMENT DEFINITION (first-order, linear):
- QUAD or QUAD4: 4 nodes (corners only), counter-clockwise order
- TRI or TRI3: 3 nodes (corners only), counter-clockwise order

ELEMENT DEFINITION (second-order, quadratic - curved edges):
- QUAD8: 8 nodes (4 corners + 4 mid-edge), counter-clockwise order
         Node order: corners 1-4, then mid-edges 5-8
- TRI6: 6 nodes (3 corners + 3 mid-edge), counter-clockwise order
        Node order: corners 1-3, then mid-edges 4-6
"""

import adsk.core
import adsk.fusion
import traceback
import json
import os
from datetime import datetime

# ============================================================================
# CONFIGURATION - Modify these values for your template
# ============================================================================

TEMPLATE_NAME = "ring_150_2D"
OUTPUT_FILENAME = "ring_template.json"

# Sketch and plane auto-detection
# The script will automatically find a sketch with this name
TEMPLATE_SKETCH_NAME = "template"
# Optional: Validate the sketch is on this plane (set to None to skip validation)
MIDPLANE_NAME = "midplane"

# Boundary node IDs (identify which nodes are at the edges after revolution)
# These must be set based on YOUR specific ring geometry after defining nodes.
# Leave as None/empty if not yet determined - you can add them to the JSON later.
#   inf_int: Node at inner edge, bottom surface (z=0, smallest x)
#   inf_ext: Node at outer edge, bottom surface (z=0, largest x)
#   sup_int: Node at inner edge, top surface (z=max, smallest x)
#   sup_ext: Node at outer edge, top surface (z=max, largest x)
BOUNDARY_NODES = {
    "inf_int": None,
    "inf_ext": None,
    "sup_int": None,
    "sup_ext": None
}

# Base height for scaling reference (mm)
# This should be the actual height (Z extent) of your ring profile.
# Set to None to auto-calculate from node coordinates, or specify manually.
BASE_HEIGHT_MM = None  # Auto-calculated from nodes if None

# Fusion 360 document units setting
# Set to "mm" if your Fusion 360 document uses millimeters
# Set to "cm" if your Fusion 360 document uses centimeters (Fusion default)
# The script will convert to mm for the output template
FUSION_UNITS = "mm"  # Change to "cm" if using Fusion's default centimeter units

# Set to True to use predefined elements instead of interactive definition
USE_PREDEFINED_ELEMENTS = False

# Predefined element connectivity (if USE_PREDEFINED_ELEMENTS is True)
# Format: {"id": int, "type": "QUAD"|"QUAD8"|"TRI"|"TRI6", "nodes": [...]}
# Define your elements here based on your specific ring geometry node IDs
#
# First-order (linear) elements - straight edges:
#   {"id": 1, "type": "QUAD", "nodes": [1, 2, 4, 3]},      # 4 corners, CCW
#   {"id": 2, "type": "TRI", "nodes": [5, 6, 7]},          # 3 corners, CCW
#
# Second-order (quadratic) elements - curved edges possible:
#   {"id": 1, "type": "QUAD8", "nodes": [1,2,3,4, 5,6,7,8]}, # 4 corners + 4 mid-edge
#   {"id": 2, "type": "TRI6", "nodes": [1,2,3, 4,5,6]},      # 3 corners + 3 mid-edge
#
# For QUAD8: nodes 1-4 are corners (CCW), nodes 5-8 are mid-edges (CCW, starting between 1-2)
# For TRI6: nodes 1-3 are corners (CCW), nodes 4-6 are mid-edges (CCW, starting between 1-2)
PREDEFINED_ELEMENTS = [
    # Add your elements here after defining nodes
]

# Surface definitions (optional, for contact surfaces in FEBio)
# Define element faces for contact surfaces if needed
# Example:
#   "RING_SURF_SUP": [{"element": 1, "face": 1}, {"element": 2, "face": 1}],
PREDEFINED_SURFACES = {
    # Add your surfaces here after defining elements
}

# ============================================================================
# MAIN SCRIPT
# ============================================================================

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)

        if not design:
            ui.messageBox("No active Fusion 360 design found.\n\nPlease open a design with a ring profile sketch.")
            return

        # Get the active sketch or let user select one
        sketch = get_target_sketch(ui, design)
        if not sketch:
            return

        # Extract nodes from sketch points
        nodes, filter_stats = extract_nodes_from_sketch(sketch)
        if not nodes:
            ui.messageBox(
                "No valid sketch points found in the selected sketch.\n\n"
                f"Points examined: {filter_stats['total_examined']}\n"
                f"Skipped (origin): {filter_stats['skipped_origin']}\n"
                f"Skipped (circle/arc centers): {filter_stats['skipped_centers']}\n"
                f"Skipped (duplicates): {filter_stats['skipped_duplicates']}\n\n"
                "Add sketch points using the Point tool to define mesh nodes."
            )
            return

        # Show filtering summary if any points were filtered
        filtered_count = filter_stats["skipped_centers"] + filter_stats["skipped_duplicates"]
        if filtered_count > 0:
            ui.messageBox(
                f"Point Extraction Summary\n\n"
                f"Total points examined: {filter_stats['total_examined']}\n"
                f"Valid nodes extracted: {len(nodes)}\n\n"
                f"Filtered out:\n"
                f"  • Circle/arc centers: {filter_stats['skipped_centers']}\n"
                f"  • Duplicate coordinates: {filter_stats['skipped_duplicates']}\n"
                f"  • Origin point: {filter_stats['skipped_origin']}",
                "Point Filtering Applied"
            )

        # Get elements (predefined or interactive)
        if USE_PREDEFINED_ELEMENTS:
            elements = PREDEFINED_ELEMENTS
            surfaces = PREDEFINED_SURFACES
        else:
            elements = define_elements_interactive(ui, nodes)
            surfaces = {}
            if not elements:
                ui.messageBox("No elements defined. Template export cancelled.")
                return

        # Validate elements reference existing nodes
        node_ids = {n["id"] for n in nodes}
        for elem in elements:
            for nid in elem["nodes"]:
                if nid not in node_ids:
                    ui.messageBox(
                        f"Element {elem['id']} references node {nid} which doesn't exist.\n\n"
                        f"Available node IDs: {sorted(node_ids)}"
                    )
                    return

        # Build the JSON structure
        template = build_template_json(nodes, elements, surfaces)

        # Save to file
        output_path = get_output_path(ui)
        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                json.dump(template, f, indent=2)
            ui.messageBox(
                f"Template exported successfully!\n\n"
                f"File: {output_path}\n"
                f"Nodes: {len(nodes)}\n"
                f"Elements: {len(elements)}\n"
                f"Surfaces: {len(surfaces)}"
            )

    except:
        if ui:
            ui.messageBox(f"Script failed:\n\n{traceback.format_exc()}")


def get_target_sketch(ui, design):
    """Find the 'template' sketch automatically, or let user select."""

    # First, try to auto-detect sketch by name
    template_sketch = None
    for comp in design.allComponents:
        for sketch in comp.sketches:
            if sketch.name.lower() == TEMPLATE_SKETCH_NAME.lower():
                # Validate it's on the correct plane if MIDPLANE_NAME is set
                if MIDPLANE_NAME:
                    ref_plane = sketch.referencePlane
                    plane_name = getattr(ref_plane, 'name', None)
                    if plane_name and plane_name.lower() != MIDPLANE_NAME.lower():
                        # Found sketch but on wrong plane, continue searching
                        continue
                template_sketch = sketch
                break
        if template_sketch:
            break

    if template_sketch:
        # Found the template sketch, confirm with user
        plane_info = ""
        if MIDPLANE_NAME:
            ref_plane = template_sketch.referencePlane
            plane_name = getattr(ref_plane, 'name', 'Unknown')
            plane_info = f"\nPlane: {plane_name}"

        result = ui.messageBox(
            f"Found sketch '{template_sketch.name}'{plane_info}\n\n"
            f"Use this sketch for template export?",
            "Template Sketch Found",
            adsk.core.MessageBoxButtonTypes.YesNoCancelButtonType
        )
        if result == adsk.core.DialogResults.DialogYes:
            return template_sketch
        elif result == adsk.core.DialogResults.DialogCancel:
            return None
        # If No, fall through to manual selection

    # Try to use active sketch if no auto-detect or user said No
    active_obj = design.activeEditObject
    if active_obj and active_obj.objectType == adsk.fusion.Sketch.classType():
        result = ui.messageBox(
            f"Use active sketch '{active_obj.name}'?",
            "Select Sketch",
            adsk.core.MessageBoxButtonTypes.YesNoCancelButtonType
        )
        if result == adsk.core.DialogResults.DialogYes:
            return active_obj
        elif result == adsk.core.DialogResults.DialogCancel:
            return None

    # Otherwise, prompt user to select manually
    ui.messageBox(
        f"Could not find sketch named '{TEMPLATE_SKETCH_NAME}'"
        + (f" on plane '{MIDPLANE_NAME}'" if MIDPLANE_NAME else "")
        + ".\n\nPlease select the sketch manually.",
        "Manual Selection Required"
    )
    selection = ui.selectEntity(
        "Select the sketch containing the ring profile points",
        "Sketches"
    )
    if selection:
        return selection.entity
    return None


def extract_nodes_from_sketch(sketch):
    """Extract all sketch points as node records.

    Points are numbered by their creation order in the sketch.
    Coordinates are converted to mm and mapped to X-Z plane (Y=0).

    Filters out:
    - Origin point (0,0,0)
    - Circle/arc center points (geometric construction points)
    - Duplicate coordinates (within tolerance)

    Returns
    -------
    tuple
        (nodes, filter_stats) where filter_stats is a dict with counts
        of skipped points by reason.
    """

    # Build set of center points to exclude (circle/arc centers)
    # These are automatically created when drawing circles/arcs and
    # are not explicitly placed points
    center_points = set()

    # Collect circle center points
    for circle in sketch.sketchCurves.sketchCircles:
        center_points.add(circle.centerSketchPoint)

    # Collect arc center points
    for arc in sketch.sketchCurves.sketchArcs:
        center_points.add(arc.centerSketchPoint)

    # Collect ellipse center points
    for ellipse in sketch.sketchCurves.sketchEllipses:
        center_points.add(ellipse.centerSketchPoint)

    # Collect elliptical arc center points
    for elliptical_arc in sketch.sketchCurves.sketchEllipticalArcs:
        center_points.add(elliptical_arc.centerSketchPoint)

    nodes = []
    node_id = 1
    seen_coords = []  # List of (x, z) tuples for duplicate detection
    tolerance = 1e-6  # Coordinate tolerance for duplicate detection (mm)

    filter_stats = {
        "skipped_origin": 0,
        "skipped_centers": 0,
        "skipped_duplicates": 0,
        "total_examined": 0,
    }

    # Get all sketch points, filtering out unwanted ones
    for point in sketch.sketchPoints:
        filter_stats["total_examined"] += 1

        # Skip the origin reference point
        if point.geometry.isEqualTo(adsk.core.Point3D.create(0, 0, 0)):
            filter_stats["skipped_origin"] += 1
            continue

        # Skip circle/arc/ellipse center points (construction geometry)
        if point in center_points:
            filter_stats["skipped_centers"] += 1
            continue

        geo = point.geometry

        # Convert coordinates based on Fusion document units
        # Map Fusion XY to template XZ (profile in XZ plane, Y=0)
        scale = 10.0 if FUSION_UNITS == "cm" else 1.0  # cm->mm or mm->mm
        x = round(geo.x * scale, 9)
        z = round(geo.y * scale, 9)  # Fusion Y -> Template Z

        # Check for duplicate coordinates (within tolerance)
        is_duplicate = False
        for seen_x, seen_z in seen_coords:
            if abs(x - seen_x) < tolerance and abs(z - seen_z) < tolerance:
                is_duplicate = True
                filter_stats["skipped_duplicates"] += 1
                break

        if is_duplicate:
            continue

        seen_coords.append((x, z))

        node = {
            "id": node_id,
            "x": x,
            "y": 0.0,
            "z": z
        }
        nodes.append(node)
        node_id += 1

    return nodes, filter_stats


def define_elements_interactive(ui, nodes):
    """Interactive element definition.

    Prompts the user to enter element connectivity manually.
    For production use with many elements, modify PREDEFINED_ELEMENTS instead.
    
    Supported element types:
    - QUAD/QUAD4: 4 nodes (linear)
    - QUAD8: 8 nodes (quadratic, curved edges)
    - TRI/TRI3: 3 nodes (linear)
    - TRI6: 6 nodes (quadratic, curved edges)
    """

    # Valid element types and their expected node counts
    ELEMENT_TYPES = {
        "QUAD": 4, "QUAD4": 4,
        "QUAD8": 8,
        "TRI": 3, "TRI3": 3,
        "TRI6": 6,
    }

    elements = []
    element_id = 1

    # Build node list display (first 15 nodes)
    node_display = []
    for n in nodes[:15]:
        node_display.append(f"  {n['id']:3d}: ({n['x']:8.4f}, {n['z']:8.4f})")
    node_list = "\n".join(node_display)
    if len(nodes) > 15:
        node_list += f"\n  ... and {len(nodes) - 15} more nodes"

    ui.messageBox(
        f"Node extraction complete: {len(nodes)} nodes found.\n\n"
        f"You will now define elements interactively.\n"
        f"Enter empty value to finish.\n\n"
        f"Format: TYPE,N1,N2,...\n"
        f"  QUAD,1,2,4,3    (4 nodes, linear)\n"
        f"  QUAD8,1,2,3,4,5,6,7,8  (8 nodes, quadratic)\n"
        f"  TRI,5,6,7       (3 nodes, linear)\n"
        f"  TRI6,1,2,3,4,5,6  (6 nodes, quadratic)",
        "Element Definition"
    )

    while True:
        # Prompt for element definition
        result = ui.inputBox(
            f"Element {element_id}\n\n"
            f"Enter: TYPE,N1,N2,...\n"
            f"(empty to finish)\n\n"
            f"Nodes:\n{node_list}",
            f"Element {element_id}",
            ""
        )

        input_value, cancelled = result

        if cancelled or not input_value.strip():
            break

        try:
            parts = [p.strip() for p in input_value.split(",")]
            elem_type = parts[0].upper()
            node_ids = [int(p) for p in parts[1:]]

            if elem_type not in ELEMENT_TYPES:
                ui.messageBox(
                    f"Unknown element type: {elem_type}\n\n"
                    f"Valid types: {', '.join(ELEMENT_TYPES.keys())}"
                )
                continue

            expected_nodes = ELEMENT_TYPES[elem_type]
            if len(node_ids) != expected_nodes:
                ui.messageBox(
                    f"{elem_type} requires {expected_nodes} nodes, got {len(node_ids)}.\n\n"
                    f"Format: {elem_type},n1,n2,...,n{expected_nodes}"
                )
                continue

            # Normalize type names (QUAD4->QUAD, TRI3->TRI, etc.)
            normalized_type = elem_type
            if elem_type == "QUAD4":
                normalized_type = "QUAD"
            elif elem_type == "TRI3":
                normalized_type = "TRI"

            elements.append({
                "id": element_id,
                "type": normalized_type,
                "nodes": node_ids
            })
            element_id += 1

        except Exception as e:
            ui.messageBox(f"Could not parse element definition:\n{e}")

    return elements


def build_template_json(nodes, elements, surfaces):
    """Construct the complete JSON template structure."""

    # Auto-calculate base height from nodes if not specified
    if BASE_HEIGHT_MM is None and nodes:
        z_coords = [n["z"] for n in nodes]
        base_height = max(z_coords) - min(z_coords)
    else:
        base_height = BASE_HEIGHT_MM or 0.0

    return {
        "name": TEMPLATE_NAME,
        "version": "1.0",
        "units": "mm",
        "profile_plane": "XZ",
        "base_height_mm": round(base_height, 6),

        "metadata": {
            "author": os.getenv("USERNAME", os.getenv("USER", "Unknown")),
            "created": datetime.now().strftime("%Y-%m-%d"),
            "source": "Fusion 360 export",
            "description": f"{TEMPLATE_NAME} - 2D cross-section profile for revolution",
            "notes": "Auto-generated from Fusion 360 sketch"
        },

        "nodes": nodes,
        "elements": elements,
        "boundary_nodes": {k: v for k, v in BOUNDARY_NODES.items() if v is not None},
        "surfaces": surfaces
    }


def get_output_path(ui):
    """Prompt user for output file location."""

    dialog = ui.createFileDialog()
    dialog.title = "Save Ring Template JSON"
    dialog.filter = "JSON Files (*.json)"
    dialog.filterIndex = 0
    dialog.initialFilename = OUTPUT_FILENAME

    if dialog.showSave() == adsk.core.DialogResults.DialogOK:
        return dialog.filename
    return None
