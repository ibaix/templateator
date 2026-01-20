"""
GMSH Mesh Generator for Ring Templates

Reads a boundary JSON file exported from Fusion 360 and generates
a structured quad mesh using GMSH.

Usage (from project root):
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json -n 10
    python scripts/mesh_with_gmsh.py boundaries/boundary.json --preview

Options:
    -o, --output    Output template JSON file (default: output_templates/template.json)
    -n, --divisions Number of mesh divisions per edge (default: 5)
    -q, --quality   Enable mesh quality optimization
    --preview       Show mesh in GMSH GUI before saving
    --element-size  Target element size in mm (alternative to divisions)
"""

import argparse
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import gmsh
except ImportError:
    print("Error: gmsh package not found.")
    print("Install with: pip install gmsh")
    sys.exit(1)


def load_boundary(filepath: str) -> dict:
    """Load boundary JSON from Fusion 360 export.
    
    Parameters
    ----------
    filepath : str
        Path to the boundary JSON file
        
    Returns
    -------
    dict
        Boundary data with curves and vertices
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        return json.load(f)


def save_template(template: dict, filepath: str):
    """Save mesh template to JSON file.
    
    Parameters
    ----------
    template : dict
        Template data with nodes and elements
    filepath : str
        Output file path
    """
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(template, f, indent=2)


def create_gmsh_geometry(boundary: dict) -> Tuple[int, List[int]]:
    """Create GMSH geometry from boundary curves.
    
    Parameters
    ----------
    boundary : dict
        Boundary data from Fusion 360
        
    Returns
    -------
    tuple
        (surface_tag, curve_tags) for the created geometry
    """
    curves = boundary["boundary"]["curves"]
    
    # Track created points to avoid duplicates
    point_tags = {}  # (x, z) -> tag
    curve_tags = []
    tolerance = 1e-6
    
    def get_or_create_point(x: float, z: float) -> int:
        """Get existing point or create new one."""
        # Round coordinates for lookup
        key = (round(x / tolerance) * tolerance, 
               round(z / tolerance) * tolerance)
        
        if key not in point_tags:
            # Create point (GMSH uses 3D, so y=0)
            tag = gmsh.model.geo.addPoint(x, z, 0)
            point_tags[key] = tag
        
        return point_tags[key]
    
    # Create curves
    for curve in curves:
        curve_type = curve["type"]
        
        if curve_type == "line":
            start = curve["start"]
            end = curve["end"]
            
            p1 = get_or_create_point(start[0], start[1])
            p2 = get_or_create_point(end[0], end[1])
            
            if p1 != p2:  # Avoid degenerate lines
                line_tag = gmsh.model.geo.addLine(p1, p2)
                curve_tags.append(line_tag)
        
        elif curve_type == "arc":
            center = curve["center"]
            radius = curve["radius"]
            start = curve["start"]
            end = curve["end"]
            
            # Create center and endpoint points
            pc = get_or_create_point(center[0], center[1])
            p1 = get_or_create_point(start[0], start[1])
            p2 = get_or_create_point(end[0], end[1])
            
            if p1 != p2:  # Avoid degenerate arcs
                # GMSH addCircleArc takes: start, center, end
                arc_tag = gmsh.model.geo.addCircleArc(p1, pc, p2)
                curve_tags.append(arc_tag)
            elif abs(curve.get("end_angle", 0) - curve.get("start_angle", 0)) > 5:
                # Full circle - need to split into two arcs
                # Create midpoint
                mid_angle = (curve.get("start_angle", 0) + math.pi) % (2 * math.pi)
                mid_x = center[0] + radius * math.cos(mid_angle)
                mid_z = center[1] + radius * math.sin(mid_angle)
                pm = get_or_create_point(mid_x, mid_z)
                
                arc1 = gmsh.model.geo.addCircleArc(p1, pc, pm)
                arc2 = gmsh.model.geo.addCircleArc(pm, pc, p1)
                curve_tags.extend([arc1, arc2])
        
        elif curve_type == "spline":
            points = curve["points"]
            if len(points) >= 2:
                # Create points for spline
                spline_points = []
                for pt in points:
                    p = get_or_create_point(pt[0], pt[1])
                    if not spline_points or p != spline_points[-1]:
                        spline_points.append(p)
                
                if len(spline_points) >= 2:
                    # Use BSpline for smooth curves
                    spline_tag = gmsh.model.geo.addBSpline(spline_points)
                    curve_tags.append(spline_tag)
    
    if not curve_tags:
        raise ValueError("No valid curves created from boundary")
    
    # Create curve loop (closed boundary)
    loop_tag = gmsh.model.geo.addCurveLoop(curve_tags)
    
    # Create plane surface
    surface_tag = gmsh.model.geo.addPlaneSurface([loop_tag])
    
    return surface_tag, curve_tags


def configure_quad_mesh(surface_tag: int, curve_tags: List[int], 
                        divisions: int = 5, element_size: Optional[float] = None):
    """Configure GMSH for quad meshing.
    
    Uses transfinite meshing for simple 3-4 sided surfaces,
    otherwise uses unstructured quad meshing with recombination.
    
    Parameters
    ----------
    surface_tag : int
        GMSH surface tag
    curve_tags : list
        List of curve tags in the boundary
    divisions : int
        Number of divisions per edge
    element_size : float, optional
        Target element size (overrides divisions if set)
    """
    # Synchronize geometry
    gmsh.model.geo.synchronize()
    
    # Set mesh options for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # Blossom
    gmsh.option.setNumber("Mesh.ElementOrder", 1)  # Linear elements
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    
    # Calculate mesh size from boundary if not specified
    if element_size is None:
        # Estimate element size from geometry bounds and divisions
        try:
            bounds = gmsh.model.getBoundingBox(-1, -1)
            diag = math.sqrt((bounds[3]-bounds[0])**2 + (bounds[4]-bounds[1])**2)
            element_size = diag / (divisions * 2)
        except:
            element_size = 0.01  # Default fallback
    
    # Set mesh size at all points
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", element_size * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", element_size * 1.5)
    
    # Check if transfinite meshing is possible (3 or 4 corners)
    num_curves = len(curve_tags)
    use_transfinite = num_curves in [3, 4]
    
    if use_transfinite:
        print(f"  Using transfinite meshing ({num_curves} curves)")
        try:
            # Set transfinite curves
            for curve_tag in curve_tags:
                try:
                    gmsh.model.mesh.setTransfiniteCurve(curve_tag, divisions + 1)
                except:
                    pass
            
            # Set transfinite surface
            gmsh.model.mesh.setTransfiniteSurface(surface_tag)
            gmsh.model.mesh.setRecombine(2, surface_tag)
        except Exception as e:
            print(f"  Transfinite failed, using unstructured: {e}")
            use_transfinite = False
    
    if not use_transfinite:
        print(f"  Using unstructured quad meshing ({num_curves} curves)")
        # Use unstructured meshing with quad recombination
        gmsh.model.mesh.setRecombine(2, surface_tag)
        
        # Set mesh size on curves for better control
        for curve_tag in curve_tags:
            try:
                gmsh.model.mesh.setTransfiniteCurve(curve_tag, divisions + 1)
            except:
                pass


def extract_mesh_data() -> Tuple[List[dict], List[dict]]:
    """Extract nodes and elements from GMSH mesh.
    
    Returns
    -------
    tuple
        (nodes, elements) lists in template format
    """
    nodes = []
    elements = []
    
    # Get all nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    
    # Create node ID mapping (GMSH uses 1-based tags)
    node_map = {}  # gmsh_tag -> our_id
    
    for i, tag in enumerate(node_tags):
        our_id = i + 1
        node_map[tag] = our_id
        
        # Coordinates are flattened [x1, y1, z1, x2, y2, z2, ...]
        x = round(node_coords[i * 3], 9)
        z = round(node_coords[i * 3 + 1], 9)  # GMSH y -> our z
        # We use y=0 for 2D templates
        
        nodes.append({
            "id": our_id,
            "x": x,
            "y": 0.0,
            "z": z
        })
    
    # Get all 2D elements (triangles and quads)
    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(2)
    
    elem_id = 1
    for i, elem_type in enumerate(element_types):
        # Get element type info
        elem_name, dim, order, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(elem_type)
        
        # Determine our element type
        if num_nodes == 4:
            our_type = "QUAD"
        elif num_nodes == 3:
            our_type = "TRI"
        elif num_nodes == 8:
            our_type = "QUAD8"
        elif num_nodes == 6:
            our_type = "TRI6"
        else:
            continue  # Skip unsupported element types
        
        # Extract each element
        node_tags_flat = element_node_tags[i]
        num_elements = len(node_tags_flat) // num_nodes
        
        for j in range(num_elements):
            elem_nodes = []
            for k in range(num_nodes):
                gmsh_node = node_tags_flat[j * num_nodes + k]
                our_node = node_map.get(gmsh_node)
                if our_node:
                    elem_nodes.append(our_node)
            
            if len(elem_nodes) == num_nodes:
                elements.append({
                    "id": elem_id,
                    "type": our_type,
                    "nodes": elem_nodes
                })
                elem_id += 1
    
    return nodes, elements


def detect_boundary_nodes(nodes: List[dict]) -> dict:
    """Detect the four corner boundary nodes.
    
    Parameters
    ----------
    nodes : list
        List of node dictionaries
        
    Returns
    -------
    dict
        Dictionary with inf_int, inf_ext, sup_int, sup_ext node IDs
    """
    if not nodes:
        return {}
    
    # Find min/max coordinates
    x_coords = [n["x"] for n in nodes]
    z_coords = [n["z"] for n in nodes]
    
    x_min, x_max = min(x_coords), max(x_coords)
    z_min, z_max = min(z_coords), max(z_coords)
    
    # Find nodes closest to corners
    def find_closest_node(target_x: float, target_z: float) -> int:
        best_id = nodes[0]["id"]
        best_dist = float('inf')
        
        for node in nodes:
            dist = math.sqrt((node["x"] - target_x)**2 + (node["z"] - target_z)**2)
            if dist < best_dist:
                best_dist = dist
                best_id = node["id"]
        
        return best_id
    
    return {
        "inf_int": find_closest_node(x_min, z_min),  # Inner edge, bottom
        "inf_ext": find_closest_node(x_max, z_min),  # Outer edge, bottom
        "sup_int": find_closest_node(x_min, z_max),  # Inner edge, top
        "sup_ext": find_closest_node(x_max, z_max),  # Outer edge, top
    }


def detect_boundary_surfaces(nodes: List[dict], elements: List[dict]) -> dict:
    """Detect boundary surfaces (top and bottom edges).
    
    Parameters
    ----------
    nodes : list
        List of node dictionaries
    elements : list
        List of element dictionaries
        
    Returns
    -------
    dict
        Dictionary with RING_SURF_SUP and RING_SURF_INF surface definitions
    """
    if not nodes or not elements:
        return {}
    
    # Find min/max Z coordinates
    z_coords = [n["z"] for n in nodes]
    z_min, z_max = min(z_coords), max(z_coords)
    z_tol = (z_max - z_min) * 0.01  # 1% tolerance
    
    # Create node lookup
    node_map = {n["id"]: n for n in nodes}
    
    # Find nodes on top and bottom boundaries
    top_nodes = {n["id"] for n in nodes if abs(n["z"] - z_max) < z_tol}
    bottom_nodes = {n["id"] for n in nodes if abs(n["z"] - z_min) < z_tol}
    
    surf_sup = []  # Top surface
    surf_inf = []  # Bottom surface
    
    for elem in elements:
        elem_nodes = elem["nodes"]
        elem_type = elem["type"]
        
        # Determine which faces are on boundaries
        # For QUAD: faces are edges 1-2, 2-3, 3-4, 4-1
        # For TRI: faces are edges 1-2, 2-3, 3-1
        
        if elem_type == "QUAD" and len(elem_nodes) >= 4:
            edges = [
                (elem_nodes[0], elem_nodes[1], 1),
                (elem_nodes[1], elem_nodes[2], 2),
                (elem_nodes[2], elem_nodes[3], 3),
                (elem_nodes[3], elem_nodes[0], 4),
            ]
        elif elem_type == "TRI" and len(elem_nodes) >= 3:
            edges = [
                (elem_nodes[0], elem_nodes[1], 1),
                (elem_nodes[1], elem_nodes[2], 2),
                (elem_nodes[2], elem_nodes[0], 3),
            ]
        else:
            continue
        
        for n1, n2, face_idx in edges:
            # Check if this edge is on top boundary
            if n1 in top_nodes and n2 in top_nodes:
                surf_sup.append({"element": elem["id"], "face": face_idx})
            # Check if this edge is on bottom boundary
            elif n1 in bottom_nodes and n2 in bottom_nodes:
                surf_inf.append({"element": elem["id"], "face": face_idx})
    
    surfaces = {}
    if surf_sup:
        surfaces["RING_SURF_SUP"] = surf_sup
    if surf_inf:
        surfaces["RING_SURF_INF"] = surf_inf
    
    return surfaces


def build_template(boundary: dict, nodes: List[dict], elements: List[dict]) -> dict:
    """Build the complete template JSON structure.
    
    Parameters
    ----------
    boundary : dict
        Original boundary data from Fusion 360
    nodes : list
        List of node dictionaries
    elements : list
        List of element dictionaries
        
    Returns
    -------
    dict
        Complete template in ring_template.json format
    """
    # Calculate base height from Z extent
    z_coords = [n["z"] for n in nodes]
    base_height = max(z_coords) - min(z_coords) if z_coords else 0.0
    
    # Detect boundary nodes and surfaces
    boundary_nodes = detect_boundary_nodes(nodes)
    surfaces = detect_boundary_surfaces(nodes, elements)
    
    # Count element types
    quad_count = sum(1 for e in elements if e["type"] == "QUAD")
    tri_count = sum(1 for e in elements if e["type"] == "TRI")
    
    return {
        "name": boundary.get("name", "ring_2D").replace("_boundary", "_2D"),
        "version": "1.0",
        "units": boundary.get("units", "mm"),
        "profile_plane": "XZ",
        "base_height_mm": round(base_height, 6),
        
        "metadata": {
            "author": os.getenv("USERNAME", os.getenv("USER", "Unknown")),
            "created": datetime.now().strftime("%Y-%m-%d"),
            "source": "Fusion 360 + GMSH",
            "description": f"{boundary.get('name', 'ring')} - 2D cross-section profile for revolution",
            "notes": f"Auto-generated mesh: {len(nodes)} nodes, {len(elements)} elements ({quad_count} quads, {tri_count} tris)"
        },
        
        "nodes": nodes,
        "elements": elements,
        "boundary_nodes": boundary_nodes,
        "surfaces": surfaces
    }


def generate_mesh(boundary_file: str, output_file: str, 
                  divisions: int = 5, element_size: Optional[float] = None,
                  quality: bool = False, preview: bool = False) -> dict:
    """Generate mesh from boundary file.
    
    Parameters
    ----------
    boundary_file : str
        Path to boundary JSON from Fusion 360
    output_file : str
        Path for output template JSON
    divisions : int
        Number of mesh divisions per edge
    element_size : float, optional
        Target element size (overrides divisions)
    quality : bool
        Enable mesh quality optimization
    preview : bool
        Show mesh in GMSH GUI before saving
        
    Returns
    -------
    dict
        Generated template data
    """
    # Load boundary
    print(f"Loading boundary: {boundary_file}")
    boundary = load_boundary(boundary_file)
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.model.add("ring_mesh")
    
    try:
        # Create geometry
        print("Creating geometry...")
        surface_tag, curve_tags = create_gmsh_geometry(boundary)
        
        # Configure meshing
        print(f"Configuring mesh (divisions={divisions})...")
        configure_quad_mesh(surface_tag, curve_tags, divisions, element_size)
        
        # Generate mesh
        print("Generating mesh...")
        gmsh.model.mesh.generate(2)
        
        # Optimize if requested
        if quality:
            print("Optimizing mesh quality...")
            gmsh.model.mesh.optimize("Laplace2D")
            gmsh.model.mesh.optimize("Relocate2D")
        
        # Preview if requested
        if preview:
            print("Opening GMSH GUI for preview...")
            gmsh.fltk.run()
        
        # Extract mesh data
        print("Extracting mesh data...")
        nodes, elements = extract_mesh_data()
        
        # Build template
        template = build_template(boundary, nodes, elements)
        
        # Save template
        print(f"Saving template: {output_file}")
        save_template(template, output_file)
        
        # Summary
        quad_count = sum(1 for e in elements if e["type"] == "QUAD")
        tri_count = sum(1 for e in elements if e["type"] == "TRI")
        
        print(f"\nMesh generation complete!")
        print(f"  Nodes: {len(nodes)}")
        print(f"  Elements: {len(elements)}")
        print(f"    - Quads: {quad_count}")
        print(f"    - Tris: {tri_count}")
        print(f"  Base height: {template['base_height_mm']:.4f} mm")
        print(f"\nOutput saved to: {output_file}")
        
        return template
        
    finally:
        gmsh.finalize()


def main():
    parser = argparse.ArgumentParser(
        description="Generate structured quad mesh from Fusion 360 boundary export",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from project root):
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json -n 10
  python scripts/mesh_with_gmsh.py boundaries/boundary.json --preview
  python scripts/mesh_with_gmsh.py boundaries/boundary.json --element-size 0.02
        """
    )
    
    parser.add_argument("boundary_file", nargs='?', default="boundaries/boundary.json",
                        help="Input boundary JSON file from Fusion 360 (default: boundaries/boundary.json)")
    parser.add_argument("-o", "--output", 
                        default="output_templates/template.json",
                        help="Output template JSON file (default: output_templates/template.json)")
    parser.add_argument("-n", "--divisions", 
                        type=int, default=5,
                        help="Number of mesh divisions per edge (default: 5)")
    parser.add_argument("--element-size", 
                        type=float, default=None,
                        help="Target element size in mm (overrides --divisions)")
    parser.add_argument("-q", "--quality", 
                        action="store_true",
                        help="Enable mesh quality optimization")
    parser.add_argument("--preview", 
                        action="store_true",
                        help="Show mesh in GMSH GUI before saving")
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.boundary_file).exists():
        print(f"Error: Boundary file not found: {args.boundary_file}")
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate mesh
    try:
        generate_mesh(
            args.boundary_file,
            args.output,
            divisions=args.divisions,
            element_size=args.element_size,
            quality=args.quality,
            preview=args.preview
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
