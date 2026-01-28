#Author-KIROSHI
#Description-Export boundary profile from 3D body intersection with midplane for GMSH meshing

# pyright: reportMissingImports=false
# type: ignore[import]

"""
Fusion 360 Add-In Script: Boundary Profile Exporter

Run the script from within Fusion 360 (Scripts and Add-Ins > Shift+S).

This script extracts the outer boundary profile by intersecting a 3D body
with a construction plane, then exports the boundary curves to JSON format
for meshing with GMSH.

WORKFLOW:
1. Ensure your design has a 3D body and a construction plane named "midplane_real"
2. Run this script from Scripts and Add-Ins (Shift+S)
3. Script automatically intersects body with plane and extracts outer profile
4. Save the boundary.json file
5. Use mesh_with_gmsh.py to generate the mesh template

COORDINATE SYSTEM:
- Fusion 360 XY plane maps to template XZ plane (Y=0)
- X = radial direction (inner to outer)
- Z = thickness direction (bottom to top)
"""

import adsk.core
import adsk.fusion
import traceback
import json
import os
import math
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

# Name for the output template
TEMPLATE_NAME = "ring_150_boundary"
OUTPUT_FILENAME = "boundary.json"

# Output directory for boundary files (relative to script location)
# The script is in: fusion360_ring_template_export/
# The boundaries folder is: ../boundaries/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BOUNDARIES_FOLDER = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "boundaries"))

# Construction plane to use for intersection
MIDPLANE_NAME = "midplane_real"

# Fusion 360 API internal units: all geometry (sketch points, radii, etc.) is
# returned in centimeters (cm), regardless of document display units.
# We export JSON in mm, so we always convert cm -> mm (scale = 10).
SCALE_CM_TO_MM = 10.0

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
            ui.messageBox(
                "No active Fusion 360 design found.\n\n"
                "Please open a design with a 3D body and construction plane."
            )
            return

        # Step 1: Find the midplane construction plane
        plane = find_plane_by_name(design, MIDPLANE_NAME)
        if not plane:
            ui.messageBox(
                f"Could not find construction plane '{MIDPLANE_NAME}'.\n\n"
                "Please create a construction plane with this name in your design."
            )
            return

        # Step 2: Get the component and bodies
        root_comp = design.rootComponent
        bodies = root_comp.bRepBodies
        
        if bodies.count == 0:
            ui.messageBox(
                "No 3D bodies found in the design.\n\n"
                "Please create a 3D body first."
            )
            return

        # If multiple bodies, let user select which one
        target_body = None
        if bodies.count == 1:
            target_body = bodies.item(0)
        else:
            # Ask user to select a body
            ui.messageBox(
                f"Found {bodies.count} bodies. Please select the body to extract profile from.",
                "Select Body"
            )
            selection = ui.selectEntity(
                "Select the body to extract profile from",
                "Bodies"
            )
            if selection:
                target_body = selection.entity
            else:
                return

        # Step 3: Create intersection sketch
        ui.messageBox(
            f"Creating intersection with plane '{MIDPLANE_NAME}'...\n\n"
            f"Body: {target_body.name}",
            "Extracting Profile"
        )

        sketch, created_new = create_intersection_sketch(root_comp, plane, target_body)
        if not sketch:
            ui.messageBox(
                "Failed to create intersection sketch.\n\n"
                "Make sure the plane intersects the body."
            )
            return

        # Step 4: Extract the outer profile
        if sketch.profiles.count == 0:
            ui.messageBox(
                "No closed profiles found from intersection.\n\n"
                "Make sure the plane passes through the body."
            )
            # Clean up if we created a new sketch
            if created_new:
                sketch.deleteMe()
            return

        outer_profile = get_outer_profile(sketch)
        if not outer_profile:
            ui.messageBox(
                "Could not determine outer profile.\n\n"
                "The intersection may not form a closed loop."
            )
            if created_new:
                sketch.deleteMe()
            return

        # Step 5: Extract curves from the outer profile
        curves = extract_boundary_curves(outer_profile, sketch)
        vertices = extract_boundary_vertices(outer_profile, sketch)

        if not curves:
            ui.messageBox(
                "No curves extracted from profile.\n\n"
                "The profile may be empty or invalid."
            )
            if created_new:
                sketch.deleteMe()
            return

        # Step 6: Build and save JSON
        boundary_data = build_boundary_json(curves, vertices)

        # Get output file path from user
        output_path = get_output_path(ui)
        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                json.dump(boundary_data, f, indent=2)
            
            # Summary
            line_count = sum(1 for c in curves if c["type"] == "line")
            arc_count = sum(1 for c in curves if c["type"] == "arc")
            spline_count = sum(1 for c in curves if c["type"] == "spline")
            
            ui.messageBox(
                f"Boundary exported successfully!\n\n"
                f"File: {output_path}\n\n"
                f"Curves: {len(curves)}\n"
                f"  • Lines: {line_count}\n"
                f"  • Arcs: {arc_count}\n"
                f"  • Splines: {spline_count}\n\n"
                f"Vertices: {len(vertices)}\n\n"
                f"Next step: Run mesh_with_gmsh.py to generate the mesh.",
                "Export Complete"
            )

        # Clean up the temporary sketch if we created it
        if created_new:
            result = ui.messageBox(
                "Delete the temporary intersection sketch?",
                "Cleanup",
                adsk.core.MessageBoxButtonTypes.YesNoButtonType
            )
            if result == adsk.core.DialogResults.DialogYes:
                sketch.deleteMe()

    except:
        if ui:
            ui.messageBox(f"Script failed:\n\n{traceback.format_exc()}")


def find_plane_by_name(design, plane_name):
    """Find a construction plane by name in the design.
    
    Searches through all components for a construction plane matching the name.
    
    Parameters
    ----------
    design : adsk.fusion.Design
        The active Fusion 360 design
    plane_name : str
        Name of the construction plane to find
        
    Returns
    -------
    adsk.fusion.ConstructionPlane or None
        The found plane, or None if not found
    """
    plane_name_lower = plane_name.lower()
    
    # Search in all components
    for comp in design.allComponents:
        for plane in comp.constructionPlanes:
            if plane.name.lower() == plane_name_lower:
                return plane
    
    # Also check root component's origin planes
    root_comp = design.rootComponent
    for plane in root_comp.originConstructionPlanes:
        if plane.name.lower() == plane_name_lower:
            return plane
    
    return None


def create_intersection_sketch(component, plane, body):
    """Create a sketch with the body's intersection on the plane.
    
    Parameters
    ----------
    component : adsk.fusion.Component
        The component to create the sketch in
    plane : adsk.fusion.ConstructionPlane
        The plane to create the sketch on
    body : adsk.fusion.BRepBody
        The body to intersect with the plane
        
    Returns
    -------
    tuple
        (sketch, created_new) where sketch is the created sketch and
        created_new indicates if a new sketch was created
    """
    # Create a new sketch on the plane
    sketches = component.sketches
    sketch = sketches.add(plane)
    sketch.name = f"_intersection_{MIDPLANE_NAME}"
    
    # Use intersectWithSketchPlane to get the cross-section
    # This projects the intersection curves onto the sketch
    try:
        sketch.intersectWithSketchPlane([body])
    except:
        # If intersectWithSketchPlane fails, try alternative approach
        # Project the body edges that lie on the plane
        try:
            sketch.project(body)
        except:
            sketch.deleteMe()
            return None, False
    
    return sketch, True


def get_outer_profile(sketch):
    """Find the outermost (largest area) profile in the sketch.
    
    For a ring cross-section, the outer profile is the largest closed loop.
    
    Parameters
    ----------
    sketch : adsk.fusion.Sketch
        The sketch containing profiles
        
    Returns
    -------
    adsk.fusion.Profile or None
        The outer profile, or None if no profiles exist
    """
    if sketch.profiles.count == 0:
        return None
    
    if sketch.profiles.count == 1:
        return sketch.profiles.item(0)
    
    # Find the profile with the largest area (outer boundary)
    largest_profile = None
    largest_area = 0
    
    for i in range(sketch.profiles.count):
        profile = sketch.profiles.item(i)
        try:
            area = profile.areaProperties().area
            if area > largest_area:
                largest_area = area
                largest_profile = profile
        except:
            # If area calculation fails, use bounding box
            pass
    
    # Fallback: if area didn't work, use bounding box
    if largest_profile is None:
        largest_bb_area = 0
        for i in range(sketch.profiles.count):
            profile = sketch.profiles.item(i)
            try:
                bb = profile.boundingBox
                bb_area = (bb.maxPoint.x - bb.minPoint.x) * (bb.maxPoint.y - bb.minPoint.y)
                if bb_area > largest_bb_area:
                    largest_bb_area = bb_area
                    largest_profile = profile
            except:
                continue
    
    return largest_profile if largest_profile else sketch.profiles.item(0)


def extract_boundary_curves(profile, sketch):
    """Extract curves from a profile's outer loop.
    
    Parameters
    ----------
    profile : adsk.fusion.Profile
        The profile to extract curves from
    sketch : adsk.fusion.Sketch
        The sketch containing the profile
        
    Returns
    -------
    list
        List of curve dictionaries with type and coordinates
    """
    curves = []
    scale = SCALE_CM_TO_MM
    
    # Get the outer loop of the profile
    outer_loop = None
    for loop in profile.profileLoops:
        if loop.isOuter:
            outer_loop = loop
            break
    
    # If no explicit outer loop, use the first one
    if outer_loop is None and profile.profileLoops.count > 0:
        outer_loop = profile.profileLoops.item(0)
    
    if outer_loop is None:
        return curves
    
    # Extract each curve in the loop
    for profile_curve in outer_loop.profileCurves:
        curve_geo = profile_curve.geometry
        curve_data = extract_curve_geometry(curve_geo, scale)
        if curve_data:
            curves.append(curve_data)
    
    return curves


def extract_curve_geometry(curve_geo, scale):
    """Extract geometry data from a curve.
    
    Parameters
    ----------
    curve_geo : adsk.core.Curve3D
        The curve geometry to extract
    scale : float
        Scale factor (cm to mm); Fusion API uses cm internally.
        
    Returns
    -------
    dict or None
        Curve data dictionary, or None if extraction failed
    """
    curve_type = curve_geo.curveType
    
    # Line segment
    if curve_type == adsk.core.Curve3DTypes.Line3DCurveType:
        line = adsk.core.Line3D.cast(curve_geo)
        if line:
            # Get start and end points
            start = line.startPoint
            end = line.endPoint
            return {
                "type": "line",
                "start": [round(start.x * scale, 9), round(start.y * scale, 9)],
                "end": [round(end.x * scale, 9), round(end.y * scale, 9)]
            }
    
    # Circular arc
    elif curve_type == adsk.core.Curve3DTypes.Arc3DCurveType:
        arc = adsk.core.Arc3D.cast(curve_geo)
        if arc:
            center = arc.center
            radius = arc.radius * scale
            start = arc.startPoint
            end = arc.endPoint
            
            # Calculate angles
            start_angle = math.atan2(start.y - center.y, start.x - center.x)
            end_angle = math.atan2(end.y - center.y, end.x - center.x)
            
            return {
                "type": "arc",
                "center": [round(center.x * scale, 9), round(center.y * scale, 9)],
                "radius": round(radius, 9),
                "start_angle": round(start_angle, 9),
                "end_angle": round(end_angle, 9),
                "start": [round(start.x * scale, 9), round(start.y * scale, 9)],
                "end": [round(end.x * scale, 9), round(end.y * scale, 9)]
            }
    
    # Circle (full circle)
    elif curve_type == adsk.core.Curve3DTypes.Circle3DCurveType:
        circle = adsk.core.Circle3D.cast(curve_geo)
        if circle:
            center = circle.center
            radius = circle.radius * scale
            return {
                "type": "arc",
                "center": [round(center.x * scale, 9), round(center.y * scale, 9)],
                "radius": round(radius, 9),
                "start_angle": 0,
                "end_angle": 2 * math.pi,
                "start": [round((center.x + circle.radius) * scale, 9), round(center.y * scale, 9)],
                "end": [round((center.x + circle.radius) * scale, 9), round(center.y * scale, 9)]
            }
    
    # Ellipse or elliptical arc
    elif curve_type in [adsk.core.Curve3DTypes.Ellipse3DCurveType, 
                        adsk.core.Curve3DTypes.EllipticalArc3DCurveType]:
        # Convert ellipse to spline points for GMSH compatibility
        evaluator = curve_geo.evaluator
        success, start_param, end_param = evaluator.getParameterExtents()
        if success:
            points = sample_curve_points(evaluator, start_param, end_param, scale, num_points=20)
            if points:
                return {
                    "type": "spline",
                    "points": points,
                    "degree": 3
                }
    
    # NURBS curve or spline
    elif curve_type == adsk.core.Curve3DTypes.NurbsCurve3DCurveType:
        nurbs = adsk.core.NurbsCurve3D.cast(curve_geo)
        if nurbs:
            evaluator = nurbs.evaluator
            success, start_param, end_param = evaluator.getParameterExtents()
            if success:
                # Sample points along the spline
                points = sample_curve_points(evaluator, start_param, end_param, scale, num_points=20)
                if points:
                    return {
                        "type": "spline",
                        "points": points,
                        "degree": min(nurbs.degree, 3)  # GMSH typically uses degree 3 max
                    }
    
    # Fallback: sample any curve type
    try:
        evaluator = curve_geo.evaluator
        success, start_param, end_param = evaluator.getParameterExtents()
        if success:
            points = sample_curve_points(evaluator, start_param, end_param, scale, num_points=20)
            if points:
                return {
                    "type": "spline",
                    "points": points,
                    "degree": 3
                }
    except:
        pass
    
    return None


def sample_curve_points(evaluator, start_param, end_param, scale, num_points=20):
    """Sample points along a curve.
    
    Parameters
    ----------
    evaluator : adsk.core.CurveEvaluator3D
        The curve evaluator
    start_param : float
        Start parameter
    end_param : float
        End parameter
    scale : float
        Scale factor (cm to mm); Fusion API uses cm internally.
    num_points : int
        Number of points to sample
        
    Returns
    -------
    list
        List of [x, z] coordinate pairs
    """
    points = []
    for i in range(num_points):
        t = start_param + (end_param - start_param) * i / (num_points - 1)
        success, point = evaluator.getPointAtParameter(t)
        if success:
            # Map Fusion XY to template XZ
            x = round(point.x * scale, 9)
            z = round(point.y * scale, 9)  # Fusion Y -> Template Z
            points.append([x, z])
    return points


def extract_boundary_vertices(profile, sketch):
    """Extract vertex positions from the profile boundary.
    
    Parameters
    ----------
    profile : adsk.fusion.Profile
        The profile to extract vertices from
    sketch : adsk.fusion.Sketch
        The sketch containing the profile
        
    Returns
    -------
    list
        List of vertex dictionaries with x, z coordinates
    """
    vertices = []
    seen_coords = set()
    scale = SCALE_CM_TO_MM
    tolerance = 1e-6
    
    # Get the outer loop of the profile
    outer_loop = None
    for loop in profile.profileLoops:
        if loop.isOuter:
            outer_loop = loop
            break
    
    if outer_loop is None and profile.profileLoops.count > 0:
        outer_loop = profile.profileLoops.item(0)
    
    if outer_loop is None:
        return vertices
    
    # Extract vertices from curve endpoints
    for profile_curve in outer_loop.profileCurves:
        sketch_entity = profile_curve.sketchEntity
        
        # Try to get start and end points
        start_point = None
        end_point = None
        
        if hasattr(sketch_entity, 'startSketchPoint'):
            start_point = sketch_entity.startSketchPoint
        if hasattr(sketch_entity, 'endSketchPoint'):
            end_point = sketch_entity.endSketchPoint
        
        for point in [start_point, end_point]:
            if point:
                geo = point.geometry
                x = round(geo.x * scale, 9)
                z = round(geo.y * scale, 9)  # Fusion Y -> Template Z
                
                # Check for duplicates
                coord_key = (round(x / tolerance) * tolerance, 
                           round(z / tolerance) * tolerance)
                if coord_key not in seen_coords:
                    seen_coords.add(coord_key)
                    vertices.append({"x": x, "z": z})
    
    return vertices


def build_boundary_json(curves, vertices):
    """Build the boundary JSON structure.
    
    Parameters
    ----------
    curves : list
        List of curve dictionaries
    vertices : list
        List of vertex dictionaries
        
    Returns
    -------
    dict
        Complete boundary JSON structure
    """
    return {
        "name": TEMPLATE_NAME,
        "version": "1.0",
        "units": "mm",
        "plane": MIDPLANE_NAME,
        
        "metadata": {
            "author": os.getenv("USERNAME", os.getenv("USER", "Unknown")),
            "created": datetime.now().strftime("%Y-%m-%d"),
            "source": "Fusion 360 body intersection",
            "description": f"{TEMPLATE_NAME} - boundary profile for GMSH meshing"
        },
        
        "boundary": {
            "curves": curves,
            "vertices": vertices
        }
    }


def get_output_path(ui):
    """Prompt user for output file location.
    
    Defaults to the boundaries/ folder in the project.
    
    Parameters
    ----------
    ui : adsk.core.UserInterface
        The Fusion 360 user interface
        
    Returns
    -------
    str or None
        Selected file path, or None if cancelled
    """
    dialog = ui.createFileDialog()
    dialog.title = "Save Boundary Profile JSON"
    dialog.filter = "JSON Files (*.json)"
    dialog.filterIndex = 0
    dialog.initialFilename = OUTPUT_FILENAME
    
    # Set initial directory to boundaries folder
    if os.path.isdir(BOUNDARIES_FOLDER):
        dialog.initialDirectory = BOUNDARIES_FOLDER
    else:
        # Create the folder if it doesn't exist
        try:
            os.makedirs(BOUNDARIES_FOLDER, exist_ok=True)
            dialog.initialDirectory = BOUNDARIES_FOLDER
        except:
            pass  # Fall back to default directory

    if dialog.showSave() == adsk.core.DialogResults.DialogOK:
        return dialog.filename
    return None
