"""
GMSH Mesh Generator for Ring Templates (100% Quad Output)

Reads a boundary JSON file exported from Fusion 360 and generates
a guaranteed 100% quad mesh using GMSH's "All-Quads Subdivision" strategy.

This approach ensures 0 triangles in output, so downstream 3D extrusion
produces pure Hex8 elements for FEBio (no Tetrahedra or Wedges).

Strategy:
    Instead of recombining triangles into quads (which can leave residual triangles),
    this script:
    1. Generates a pure triangular mesh using Frontal 2D algorithm
    2. Subdivides each triangle into 3 quads mathematically
    
    The subdivision is topologically guaranteed to produce only quads.

Usage (from project root):
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json -n 20
    python scripts/mesh_with_gmsh.py boundaries/boundary.json --preview
    
    # High quality mesh (recommended for FEM):
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -q
    
    # Maximum quality for critical simulations:
    python scripts/mesh_with_gmsh.py boundaries/boundary.json -n 30 -q --optimize-iter 30

Options:
    -o, --output      Output template JSON file (default: output_templates/template.json)
    -n, --divisions   Number of mesh divisions per edge (default: 10)
    -q, --quality     Enable aggressive quality optimization (multiple smoothing passes)
    --optimize-iter   Number of optimization passes when -q is used (default: 40)
    --preview         Show mesh in GMSH GUI before saving
    --element-size    Target element size in mm (alternative to divisions)

Quality Metrics:
    The script reports mesh quality including:
    - Angle range: internal angles of quad elements (ideal: 90°, acceptable: 60°-120°)
    - Quality score: 0-1 metric based on angles and aspect ratio
    - Poor elements: quads with angles >150° or <30° (problematic for FEM)
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


def compute_curve_length(curve_tag: int, num_samples: int = 50) -> float:
    """Compute the length of a curve by sampling points along it.
    
    Works with the built-in geo kernel (unlike getMass which requires occ).
    
    Parameters
    ----------
    curve_tag : int
        GMSH curve tag
    num_samples : int
        Number of sample points for length estimation (more = more accurate)
        
    Returns
    -------
    float
        Approximate curve length
    """
    try:
        # Get parametric bounds of the curve
        # Returns (min_array, max_array) - for 1D curves, each has 1 element
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min = bounds[0][0]  # First element of min_array
        t_max = bounds[1][0]  # First element of max_array
        
        # Sample points along the curve
        total_length = 0.0
        prev_point = None
        
        for i in range(num_samples + 1):
            t = t_min + (t_max - t_min) * i / num_samples
            # getValue returns [x, y, z] for the point at parameter t
            point = gmsh.model.getValue(1, curve_tag, [t])
            
            if prev_point is not None:
                # Calculate distance between consecutive points
                dx = point[0] - prev_point[0]
                dy = point[1] - prev_point[1]
                dz = point[2] - prev_point[2]
                total_length += math.sqrt(dx*dx + dy*dy + dz*dz)
            
            prev_point = point
        
        return total_length
    except Exception as e:
        raise RuntimeError(f"Failed to compute curve length: {e}")


def configure_quad_mesh(surface_tag: int, curve_tags: List[int], 
                        divisions: int = 5, element_size: Optional[float] = None,
                        num_threads: int = 0, force_adaptive: bool = False):
    """Configure GMSH for 100% quad meshing using subdivision strategy.
    
    Instead of recombining triangles into quads (which can leave residual triangles),
    this uses the "All-Quads Subdivision" approach:
    1. Generate a pure triangular mesh using Delaunay algorithm
    2. Subdivide each triangle into 3 quads mathematically
    
    This guarantees 100% quad output with no residual triangles, ensuring
    downstream 3D extrusion produces pure Hex8 elements.
    
    For non-transfinite surfaces, uses GMSH's native curvature-adaptive meshing
    which automatically refines curved regions and coarsens straight edges,
    bounded by CharacteristicLengthMin/Max.
    
    Parameters
    ----------
    surface_tag : int
        GMSH surface tag
    curve_tags : list
        List of curve tags in the boundary
    divisions : int
        Number of divisions per edge (used for transfinite meshing only)
    element_size : float, optional
        Target element size (overrides divisions if set)
    num_threads : int
        Number of CPU threads (0 = auto-detect)
    force_adaptive : bool
        If True, force adaptive discretization even for 3 or 4-sided shapes.
        Useful when edges have very different lengths (e.g. ring profiles).
    """
    # Synchronize geometry (required before getMass and other queries)
    gmsh.model.geo.synchronize()
    
    # === PARALLELIZATION SETTINGS ===
    # Enable multi-threading for meshing operations
    if num_threads > 0:
        gmsh.option.setNumber("General.NumThreads", num_threads)
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", num_threads)
    
    # === MESH ALGORITHM SETTINGS (All-Quads Subdivision Strategy) ===
    # Step 1: Generate triangular mesh using Delaunay algorithm
    # Delaunay (5) handles non-symmetric shapes better than Frontal 2D (6)
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay
    
    # Step 2: Disable recombination (we use subdivision instead)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    
    # Step 3: Enable All-Quads subdivision
    # SubdivisionAlgorithm = 1 splits each triangle into 3 quads
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # All Quads
    
    gmsh.option.setNumber("Mesh.ElementOrder", 1)  # Linear elements
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    
    # === QUALITY SETTINGS ===
    # Quality type for optimization
    gmsh.option.setNumber("Mesh.QualityType", 2)  # SICN (signed inverse condition number) - best for FEM
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.3)  # Optimize elements below this quality
    
    # === SMOOTHING SETTINGS ===
    # Smoothing is beneficial for subdivided meshes (Lloyd optimization)
    gmsh.option.setNumber("Mesh.Smoothing", 1000)  # Number of smoothing steps
    gmsh.option.setNumber("Mesh.SmoothRatio", 2.0)  # Max ratio for smoothing
    gmsh.option.setNumber("Mesh.SmoothNormals", 1)  # Smooth normals for better quality
    
    # === CURVATURE-ADAPTIVE MESHING ===
    # Enable automatic element sizing based on local curvature (radius R):
    # h ~= 2 * pi * R / N_nodes
    # This refines curved regions (fillets) and coarsens straight edges automatically
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
    
    # === DENSITY ADJUSTMENT FOR SUBDIVISION ===
    # Since subdivision splits edges in half (and triangles into 3 quads),
    # we need to use a coarser pre-subdivision mesh to achieve the user's
    # expected final resolution. Factor of 2.0 compensates for edge splitting.
    subdivision_factor = 2.0
    
    adjusted_divisions = max(2, int(divisions / subdivision_factor))
    
    # Calculate target element size for adaptive discretization
    if element_size is not None:
        # User provided element_size - account for subdivision
        target_element_size = element_size * subdivision_factor
    else:
        # Estimate element size from geometry bounds and divisions
        try:
            bounds = gmsh.model.getBoundingBox(-1, -1)
            diag = math.sqrt((bounds[3]-bounds[0])**2 + (bounds[4]-bounds[1])**2)
            target_element_size = diag / (adjusted_divisions * 2)
        except:
            target_element_size = 0.01  # Default fallback
    
    # Set mesh size constraints
    # Min: ensures resolution on sharp fillets/curves (curvature field can go this small)
    # Max: allows growth on straight edges and interior (up to 5x larger)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", target_element_size * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", target_element_size * 5.0)
    
    # Enable extending mesh size from boundary into interior
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
    
    # Set mesh size at all geometry points to provide baseline sizing
    # This ensures reasonable element sizes even for straight-line-only geometries
    # where curvature-adaptive has no effect
    try:
        # Get all points in the model
        points = gmsh.model.getEntities(0)  # dim=0 for points
        point_tags = [(0, p[1]) for p in points]  # List of (dim, tag) tuples
        if point_tags:
            gmsh.model.mesh.setSize(point_tags, target_element_size)
            print(f"  Set mesh size {target_element_size:.6f} at {len(point_tags)} geometry points")
    except Exception as e:
        print(f"  Warning: Could not set point mesh sizes: {e}")
    
    # Check if transfinite meshing is possible (3 or 4 corners)
    # Can be overridden by force_adaptive for shapes with varying edge lengths
    num_curves = len(curve_tags)
    use_transfinite = (num_curves in [3, 4]) and not force_adaptive
    
    if force_adaptive and num_curves in [3, 4]:
        print(f"  Forcing adaptive discretization for {num_curves}-sided shape")
    
    if use_transfinite:
        # For transfinite surfaces, use fixed divisions to maintain grid structure
        # (opposite sides must have matching point counts)
        print(f"  Using transfinite meshing with subdivision ({num_curves} curves)")
        print(f"  Adjusted divisions: {divisions} -> {adjusted_divisions} (pre-subdivision)")
        try:
            # Set transfinite curves with adjusted divisions
            for curve_tag in curve_tags:
                try:
                    gmsh.model.mesh.setTransfiniteCurve(curve_tag, adjusted_divisions + 1)
                except:
                    pass
            
            # Set transfinite surface (no recombination needed - subdivision handles quads)
            gmsh.model.mesh.setTransfiniteSurface(surface_tag)
        except Exception as e:
            print(f"  Transfinite failed, using unstructured: {e}")
            use_transfinite = False
    
    if not use_transfinite:
        # === ADAPTIVE UNSTRUCTURED MESHING WITH BOUNDARY REFINEMENT ===
        # For contact simulations: fine elements on boundary, coarser in interior
        # Uses Distance + Threshold fields to create size gradient
        print(f"  Using adaptive unstructured meshing ({num_curves} curves)")
        print(f"  Boundary element size: {target_element_size:.6f}")
        print(f"  Interior element size: up to {target_element_size * 4.0:.6f}")
        
        try:
            # Create Distance field from all boundary curves
            # This computes distance from any point to the nearest boundary curve
            dist_field = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist_field, "CurvesList", curve_tags)
            gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 100)
            
            # Create Threshold field for size gradient
            # - Near boundary (distance < DistMin): use SizeMin (fine for contact)
            # - Far from boundary (distance > DistMax): use SizeMax (coarse interior)
            # - In between: linear interpolation
            thresh_field = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
            gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", target_element_size)
            gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", target_element_size * 8.0)
            
            # Estimate transition distance based on geometry size
            # DistMin: start growing after this distance from boundary
            # DistMax: reach maximum size at this distance
            try:
                bounds = gmsh.model.getBoundingBox(-1, -1)
                diag = math.sqrt((bounds[3]-bounds[0])**2 + (bounds[4]-bounds[1])**2)
                dist_min = diag * 0.0  # Start growing at 5% of diagonal
                dist_max = diag * 0.2   # Reach max at 30% of diagonal
            except:
                dist_min = target_element_size * 0.0
                dist_max = target_element_size * 0.2
            
            gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", dist_min)
            gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", dist_max)
            
            # Combine with curvature field using Min (take the smaller of both)
            # This ensures curved boundaries still get refined
            min_field = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [thresh_field])
            
            # Set as background mesh
            gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
            
            print(f"  Distance field: boundary refinement with gradient to interior")
            print(f"  Transition: {dist_min:.6f} - {dist_max:.6f} from boundary")
            
        except Exception as e:
            print(f"  Warning: Could not set up distance field: {e}")
            print(f"  Falling back to uniform sizing")


def compute_element_quality(nodes: List[dict], elements: List[dict]) -> dict:
    """Compute quality metrics for all elements.
    
    Calculates:
    - Internal angles (min, max per element)
    - Aspect ratio
    - Jacobian-based quality metric
    
    Parameters
    ----------
    nodes : list
        List of node dictionaries
    elements : list
        List of element dictionaries
        
    Returns
    -------
    dict
        Quality statistics and list of poor quality elements
    """
    node_map = {n["id"]: (n["x"], n["z"]) for n in nodes}
    
    all_min_angles = []
    all_max_angles = []
    all_jacobians = []
    poor_elements = []  # Elements with quality issues
    
    for elem in elements:
        if elem["type"] != "QUAD":
            continue
            
        node_ids = elem["nodes"]
        if len(node_ids) < 4:
            continue
            
        # Get coordinates
        coords = [node_map.get(nid, (0, 0)) for nid in node_ids]
        
        # Calculate internal angles at each corner
        angles = []
        for i in range(4):
            p_prev = coords[(i - 1) % 4]
            p_curr = coords[i]
            p_next = coords[(i + 1) % 4]
            
            # Vectors from current to prev and next
            v1 = (p_prev[0] - p_curr[0], p_prev[1] - p_curr[1])
            v2 = (p_next[0] - p_curr[0], p_next[1] - p_curr[1])
            
            # Magnitudes
            mag1 = math.sqrt(v1[0]**2 + v1[1]**2)
            mag2 = math.sqrt(v2[0]**2 + v2[1]**2)
            
            if mag1 < 1e-10 or mag2 < 1e-10:
                angles.append(180.0)  # Degenerate
                continue
            
            # Dot product and angle
            dot = v1[0] * v2[0] + v1[1] * v2[1]
            cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
            angle_deg = math.degrees(math.acos(cos_angle))
            angles.append(angle_deg)
        
        min_angle = min(angles)
        max_angle = max(angles)
        all_min_angles.append(min_angle)
        all_max_angles.append(max_angle)
        
        # Simplified Jacobian quality: ratio of min/max edge lengths and angle deviation
        edges = []
        for i in range(4):
            p1 = coords[i]
            p2 = coords[(i + 1) % 4]
            edge_len = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
            edges.append(edge_len)
        
        if max(edges) > 1e-10:
            aspect_ratio = max(edges) / max(min(edges), 1e-10)
        else:
            aspect_ratio = 1.0
        
        # Quality metric: 1.0 = ideal (all 90° angles), 0.0 = worst
        # Penalize angles far from 90° and high aspect ratios
        angle_quality = 1.0 - (max(abs(max_angle - 90), abs(90 - min_angle)) / 90.0)
        aspect_quality = min(1.0, 2.0 / aspect_ratio)
        jacobian_quality = angle_quality * aspect_quality
        all_jacobians.append(jacobian_quality)
        
        # Flag poor quality elements
        if max_angle > 150 or min_angle < 30 or aspect_ratio > 5:
            poor_elements.append({
                "id": elem["id"],
                "min_angle": round(min_angle, 1),
                "max_angle": round(max_angle, 1),
                "aspect_ratio": round(aspect_ratio, 2),
                "quality": round(jacobian_quality, 3)
            })
    
    # Compute statistics
    stats = {}
    if all_min_angles:
        stats["angle_min"] = round(min(all_min_angles), 1)
        stats["angle_max"] = round(max(all_max_angles), 1)
        stats["angle_avg_min"] = round(sum(all_min_angles) / len(all_min_angles), 1)
        stats["angle_avg_max"] = round(sum(all_max_angles) / len(all_max_angles), 1)
    if all_jacobians:
        stats["quality_min"] = round(min(all_jacobians), 3)
        stats["quality_max"] = round(max(all_jacobians), 3)
        stats["quality_avg"] = round(sum(all_jacobians) / len(all_jacobians), 3)
    
    return {
        "statistics": stats,
        "poor_elements": poor_elements,
        "num_poor": len(poor_elements)
    }


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
                  quality: bool = False, preview: bool = False,
                  optimize_iterations: int = 10,
                  num_threads: int = 0, force_adaptive: bool = False) -> dict:
    """Generate 100% quad mesh from boundary file using subdivision strategy.
    
    Uses the "All-Quads Subdivision" approach to guarantee pure quad output:
    1. Generates a triangular mesh
    2. Subdivides each triangle into 3 quads
    
    This ensures 0 triangles in output, so downstream 3D extrusion produces
    pure Hex8 elements for FEBio.
    
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
    optimize_iterations : int
        Number of optimization passes when quality=True
    num_threads : int
        Number of CPU threads (0 = auto-detect all cores)
    force_adaptive : bool
        Force adaptive curve discretization even for 3 or 4-sided shapes.
        Recommended for ring profiles with varying edge lengths.
        
    Returns
    -------
    dict
        Generated template data
    """
    import time
    start_time = time.perf_counter()
    
    # Load boundary
    print(f"Loading boundary: {boundary_file}", flush=True)
    boundary = load_boundary(boundary_file)
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.model.add("ring_mesh")
    
    try:
        # Create geometry
        print("Creating geometry...")
        surface_tag, curve_tags = create_gmsh_geometry(boundary)
        
        # Configure meshing
        print(f"Configuring mesh (divisions={divisions}, subdivision strategy for 100% quads)...")
        configure_quad_mesh(surface_tag, curve_tags, divisions, element_size, num_threads, force_adaptive)
        
        # Generate mesh (triangular first, then subdivide to quads)
        print("Generating mesh...", flush=True)
        mesh_start = time.perf_counter()
        
        # Step 1: Generate triangular mesh
        gmsh.model.mesh.generate(2)
        tri_time = time.perf_counter() - mesh_start
        print(f"  Triangular mesh: {tri_time:.3f}s", flush=True)
        
        # Step 2: Apply subdivision to convert triangles to quads
        # SubdivisionAlgorithm=1 splits each triangle into 3 quads
        print("  Subdividing triangles to quads...", flush=True)
        subdiv_start = time.perf_counter()
        gmsh.model.mesh.refine()
        subdiv_time = time.perf_counter() - subdiv_start
        print(f"  Subdivision: {subdiv_time:.3f}s", flush=True)
        
        mesh_time = time.perf_counter() - mesh_start
        print(f"  Total mesh generation: {mesh_time:.3f}s", flush=True)
        
        # Optimize if requested
        if quality:
            print(f"Optimizing mesh quality ({optimize_iterations} iterations)...", flush=True)
            opt_start = time.perf_counter()
            
            for i in range(optimize_iterations):
                # Laplace smoothing - moves nodes to improve element quality
                gmsh.model.mesh.optimize("Laplace2D")
                # Relocate - more aggressive node relocation  
                gmsh.model.mesh.optimize("Relocate2D")
            
            opt_time = time.perf_counter() - opt_start
            print(f"  Optimization: {opt_time:.3f}s ({optimize_iterations} iterations)", flush=True)
        
        # Preview if requested
        if preview:
            print("Opening GMSH GUI for preview...")
            gmsh.fltk.run()
        
        # Extract mesh data
        print("Extracting mesh data...", flush=True)
        try:
            nodes, elements = extract_mesh_data()
        except Exception as e:
            print(f"Error extracting mesh data: {e}", flush=True)
            raise
        
        # Compute and report quality metrics
        print("\nAnalyzing mesh quality...", flush=True)
        try:
            quality_data = compute_element_quality(nodes, elements)
            stats = quality_data["statistics"]
        except Exception as e:
            print(f"Error computing quality: {e}", flush=True)
            quality_data = {"statistics": {}, "poor_elements": [], "num_poor": 0}
            stats = {}
        
        if stats:
            print(f"  Angle range: {stats.get('angle_min', 'N/A')}° - {stats.get('angle_max', 'N/A')}°")
            print(f"  Avg angles: min={stats.get('angle_avg_min', 'N/A')}°, max={stats.get('angle_avg_max', 'N/A')}°")
            print(f"  Quality: min={stats.get('quality_min', 'N/A')}, avg={stats.get('quality_avg', 'N/A')}, max={stats.get('quality_max', 'N/A')}")
            
            if quality_data["num_poor"] > 0:
                print(f"\n  WARNING: {quality_data['num_poor']} elements with poor quality detected!")
                print(f"  (angles >150° or <30°, or aspect ratio >5)")
                # Show up to 5 worst elements
                worst = sorted(quality_data["poor_elements"], key=lambda x: x["quality"])[:5]
                for elem in worst:
                    print(f"    Element {elem['id']}: angles={elem['min_angle']}°-{elem['max_angle']}°, "
                          f"aspect={elem['aspect_ratio']}, quality={elem['quality']}")
        
        # Build template
        template = build_template(boundary, nodes, elements)
        
        # Add quality metrics to template metadata
        template["metadata"]["mesh_quality"] = stats
        if quality_data["num_poor"] > 0:
            template["metadata"]["poor_quality_elements"] = quality_data["num_poor"]
        
        # Save template
        print(f"\nSaving template: {output_file}", flush=True)
        save_template(template, output_file)
        print("Template saved successfully.", flush=True)
        
        # Summary
        total_time = time.perf_counter() - start_time
        quad_count = sum(1 for e in elements if e["type"] == "QUAD")
        tri_count = sum(1 for e in elements if e["type"] == "TRI")
        
        print(f"\n{'='*50}")
        print(f"MESH GENERATION COMPLETE")
        print(f"{'='*50}")
        print(f"  Nodes: {len(nodes)}")
        print(f"  Elements: {len(elements)}")
        print(f"    - Quads: {quad_count}")
        print(f"    - Tris: {tri_count}")
        print(f"  Base height: {template['base_height_mm']:.4f} mm")
        print(f"  Total time: {total_time:.3f}s")
        print(f"\nOutput saved to: {output_file}")
        
        # Recommendations for poor quality meshes
        if quality_data["num_poor"] > 0:
            print("\n--- RECOMMENDATIONS ---")
            print("To improve mesh quality, try:")
            print("  1. Use -q flag for quality optimization")
            print("  2. Increase --optimize-iter (e.g., 50)")
            print("  3. Increase divisions (-n) for finer mesh")
        
        return template
        
    finally:
        gmsh.finalize()


def main():
    parser = argparse.ArgumentParser(
        description="Generate 100%% quad mesh from Fusion 360 boundary using subdivision strategy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from project root):
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json -n 20
  python scripts/mesh_with_gmsh.py boundaries/boundary.json --preview
  
  # High quality mesh with optimization (recommended):
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -q
  
  # Maximum quality for critical simulations:
  python scripts/mesh_with_gmsh.py boundaries/boundary.json -n 30 -q --optimize-iter 30

Subdivision Strategy:
  This script guarantees 100%% quad output by using GMSH's "All-Quads Subdivision"
  approach, which subdivides each triangle into 3 quads. This ensures 0 residual
  triangles, so downstream 3D extrusion produces pure Hex8 elements.
        """
    )
    
    parser.add_argument("boundary_file", nargs='?', default="boundaries/boundary.json",
                        help="Input boundary JSON file from Fusion 360 (default: boundaries/boundary.json)")
    parser.add_argument("-o", "--output", 
                        default="output_templates/template.json",
                        help="Output template JSON file (default: output_templates/template.json)")
    parser.add_argument("-n", "--divisions", 
                        type=int, default=4,
                        help="Number of mesh divisions per edge (default: 10)")
    parser.add_argument("--element-size", 
                        type=float, default=None,
                        help="Target element size in mm (overrides --divisions)")
    parser.add_argument("-q", "--quality", 
                        action="store_true",
                        help="Enable aggressive quality optimization (multiple Laplace + Relocate passes)")
    parser.add_argument("--optimize-iter",
                        type=int, default=80,
                        help="Number of optimization iterations when -q is used (default: 40)")
    parser.add_argument("--threads",
                        type=int, default=0,
                        help="Number of CPU threads (default: 0 = use GMSH default)")
    parser.add_argument("--preview", 
                        action="store_true",
                        help="Show mesh in GMSH GUI before saving")
    parser.add_argument("--force-adaptive", 
                        action="store_true",
                        help="Force adaptive curve discretization even for 3 or 4-sided shapes. "
                             "Recommended for ring profiles with varying edge lengths.")
    
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
            preview=args.preview,
            optimize_iterations=args.optimize_iter,
            num_threads=args.threads,
            force_adaptive=args.force_adaptive
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
