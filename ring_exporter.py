"""Exportadores de mallas a formatos Abaqus (.inp) y FEBio (.feb)."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .ring_revolver import RevolvedMesh

# Try to import pyfebio for enhanced FEBio export
try:
    import pyfebio
    PYFEBIO_AVAILABLE = True
except ImportError:
    pyfebio = None
    PYFEBIO_AVAILABLE = False


def _format_float(value: float, precision: int) -> str:
    return f"{value:.{precision}f}"


def _ensure_parent(path: Path) -> None:
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)


def write_inp(
    mesh: RevolvedMesh,
    output_path: str | Path,
    *,
    element_type: str = "C3D8",
    precision: int = 6,
    surfaces: Optional[Dict[str, Tuple[Tuple[int, int], ...]]] = None,
) -> None:
    """Escribe la malla revolvida en formato Abaqus .inp."""

    path = Path(output_path)
    _ensure_parent(path)

    with path.open("w", encoding="utf-8") as handler:
        handler.write("*Heading\n")
        handler.write("Generated mesh by revolution script\n")
        handler.write("*Node\n")

        for node_id in sorted(mesh.nodes):
            node = mesh.nodes[node_id]
            x = _format_float(node.x, precision)
            y = _format_float(node.y, precision)
            z = _format_float(node.z, precision)
            handler.write(f"{node.id}, {x}, {y}, {z}\n")

        if mesh.hex_elements:
            # Group elements by type (C3D8, C3D20, etc.)
            hex_by_type: Dict[str, List[int]] = {}
            for elem_id in mesh.hex_elements:
                elem_type = mesh.hex_elements[elem_id].element_type
                hex_by_type.setdefault(elem_type, []).append(elem_id)

            for elem_type, elem_ids in hex_by_type.items():
                handler.write(f"*Element, type={elem_type}\n")
                for element_id in sorted(elem_ids):
                    element = mesh.hex_elements[element_id]
                    conn = ", ".join(str(nid) for nid in element.connectivity)
                    handler.write(f"{element.id}, {conn}\n")

        if mesh.prism_elements:
            # Group elements by type (C3D6, C3D15, etc.)
            prism_by_type: Dict[str, List[int]] = {}
            for elem_id in mesh.prism_elements:
                elem_type = mesh.prism_elements[elem_id].element_type
                prism_by_type.setdefault(elem_type, []).append(elem_id)

            for elem_type, elem_ids in prism_by_type.items():
                handler.write(f"*Element, type={elem_type}\n")
                for element_id in sorted(elem_ids):
                    element = mesh.prism_elements[element_id]
                    conn = ", ".join(str(nid) for nid in element.connectivity)
                    handler.write(f"{element.id}, {conn}\n")

        # Write surface definitions if provided
        if surfaces:
            for surf_name, facets in surfaces.items():
                handler.write(f"*Surface, name={surf_name}, type=ELEMENT\n")
                for elem_id, face_idx in facets:
                    handler.write(f"{elem_id}, S{face_idx}\n")

        handler.write("*End Part\n")


def write_feb(
    mesh: RevolvedMesh,
    output_path: str | Path,
    *,
    precision: int = 6,
    surfaces: Optional[Dict[str, Tuple[Tuple[int, int], ...]]] = None,
    material_name: str = "RingMaterial",
) -> None:
    """Escribe la malla revolvida en formato FEBio (.feb).
    
    Uses pyfebio if available for better compatibility, otherwise
    falls back to direct XML generation.
    """
    path = Path(output_path)
    _ensure_parent(path)

    if PYFEBIO_AVAILABLE:
        _write_feb_pyfebio(mesh, path, precision, surfaces, material_name)
    else:
        _write_feb_xml(mesh, path, precision, surfaces)


def _write_feb_pyfebio(
    mesh: RevolvedMesh,
    path: Path,
    precision: int,
    surfaces: Optional[Dict[str, Tuple[Tuple[int, int], ...]]] = None,
    material_name: str = "RingMaterial",
) -> None:
    """Write FEBio file using pyfebio library."""
    import numpy as np
    
    model = pyfebio.model.Model()
    
    # Add nodes
    nodes_section = pyfebio.mesh.Nodes(name="RingNodes")
    for node_id in sorted(mesh.nodes):
        node = mesh.nodes[node_id]
        coord_str = f"{node.x:.{precision}f},{node.y:.{precision}f},{node.z:.{precision}f}"
        nodes_section.add_node(pyfebio.mesh.Node(id=int(node.id), text=coord_str))
    model.mesh_.nodes.append(nodes_section)
    
    # Add material (Neo-Hookean default)
    material = pyfebio.material.NeoHookean(
        id=1,
        name=material_name,
        E=pyfebio.material.MaterialParameter(text=1.0),
        v=pyfebio.material.MaterialParameter(text=0.45)
    )
    model.material_.add_material(material)
    
    # Map element types to FEBio types
    FEBIO_HEX_TYPES = {"C3D8": "hex8", "C3D20": "hex20"}
    FEBIO_PRISM_TYPES = {"C3D6": "penta6", "C3D15": "penta15"}

    # Add hex elements (grouped by type)
    if mesh.hex_elements:
        hex_by_type: Dict[str, List[int]] = {}
        for elem_id in mesh.hex_elements:
            elem_type = mesh.hex_elements[elem_id].element_type
            hex_by_type.setdefault(elem_type, []).append(elem_id)

        for elem_type, elem_ids in hex_by_type.items():
            febio_type = FEBIO_HEX_TYPES.get(elem_type, "hex8")
            section_name = f"RingHex_{febio_type}"
            hex_section = pyfebio.mesh.Elements(name=section_name, type=febio_type)

            for element_id in sorted(elem_ids):
                element = mesh.hex_elements[element_id]
                conn_str = ",".join(str(nid) for nid in element.connectivity)
                # Use generic element addition for both hex8 and hex20
                if febio_type == "hex8" and hasattr(pyfebio.mesh, 'Hex8Element'):
                    hex_section.add_element(pyfebio.mesh.Hex8Element(id=int(element.id), text=conn_str))
                elif febio_type == "hex20" and hasattr(pyfebio.mesh, 'Hex20Element'):
                    hex_section.add_element(pyfebio.mesh.Hex20Element(id=int(element.id), text=conn_str))
                else:
                    # Fallback: try generic element
                    hex_section.add_element(pyfebio.mesh.Element(id=int(element.id), text=conn_str))

            model.mesh_.elements.append(hex_section)

            # Create solid domain
            domain = pyfebio.meshdomains.SolidDomain(name=section_name, mat=material_name)
            model.meshdomains_.add_solid_domain(domain)

    # Add prism elements (grouped by type)
    if mesh.prism_elements:
        prism_by_type: Dict[str, List[int]] = {}
        for elem_id in mesh.prism_elements:
            elem_type = mesh.prism_elements[elem_id].element_type
            prism_by_type.setdefault(elem_type, []).append(elem_id)

        for elem_type, elem_ids in prism_by_type.items():
            febio_type = FEBIO_PRISM_TYPES.get(elem_type, "penta6")
            section_name = f"RingPrism_{febio_type}"
            prism_section = pyfebio.mesh.Elements(name=section_name, type=febio_type)

            for element_id in sorted(elem_ids):
                element = mesh.prism_elements[element_id]
                conn_str = ",".join(str(nid) for nid in element.connectivity)
                # Use generic element addition for both penta6 and penta15
                if febio_type == "penta6" and hasattr(pyfebio.mesh, 'Penta6Element'):
                    prism_section.add_element(pyfebio.mesh.Penta6Element(id=int(element.id), text=conn_str))
                elif febio_type == "penta15" and hasattr(pyfebio.mesh, 'Penta15Element'):
                    prism_section.add_element(pyfebio.mesh.Penta15Element(id=int(element.id), text=conn_str))
                else:
                    # Fallback: try generic element
                    prism_section.add_element(pyfebio.mesh.Element(id=int(element.id), text=conn_str))

            model.mesh_.elements.append(prism_section)

            # Create solid domain
            domain = pyfebio.meshdomains.SolidDomain(name=section_name, mat=material_name)
            model.meshdomains_.add_solid_domain(domain)
    
    # Add surfaces if provided
    if surfaces:
        for surf_name, facets in surfaces.items():
            try:
                surface = pyfebio.mesh.Surface(name=surf_name)
                for elem_id, face_idx in facets:
                    if hasattr(surface, 'add_facet'):
                        surface.add_facet(elem_id, face_idx)
                    elif hasattr(surface, 'facets'):
                        surface.facets.append((elem_id, face_idx))
                model.mesh_.surfaces.append(surface)
            except (AttributeError, TypeError):
                # pyfebio surface API may be limited
                pass
    
    model.save(str(path))
    print(f"  Ring mesh exported to FEBio: {path}")


def _write_feb_xml(
    mesh: RevolvedMesh,
    path: Path,
    precision: int,
    surfaces: Optional[Dict[str, Tuple[Tuple[int, int], ...]]] = None,
) -> None:
    """Write FEBio file using direct XML generation (fallback)."""
    with path.open("w", encoding="utf-8") as handler:
        handler.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        handler.write("<febio_spec version=\"4.0\">\n")
        
        # Module
        handler.write("  <Module type=\"solid\"/>\n")
        
        # Material
        handler.write("  <Material>\n")
        handler.write("    <material id=\"1\" name=\"RingMaterial\" type=\"neo-Hookean\">\n")
        handler.write("      <E>1.0</E>\n")
        handler.write("      <v>0.45</v>\n")
        handler.write("    </material>\n")
        handler.write("  </Material>\n")
        
        # Mesh
        handler.write("  <Mesh>\n")
        handler.write("    <Nodes name=\"RingNodes\">\n")

        for node_id in sorted(mesh.nodes):
            node = mesh.nodes[node_id]
            x = _format_float(node.x, precision)
            y = _format_float(node.y, precision)
            z = _format_float(node.z, precision)
            handler.write(f"      <node id=\"{node.id}\">{x},{y},{z}</node>\n")

        handler.write("    </Nodes>\n")

        # Map element types to FEBio types
        FEBIO_HEX_TYPES = {"C3D8": "hex8", "C3D20": "hex20"}
        FEBIO_PRISM_TYPES = {"C3D6": "penta6", "C3D15": "penta15"}

        if mesh.hex_elements:
            # Group by type
            hex_by_type: Dict[str, List[int]] = {}
            for elem_id in mesh.hex_elements:
                elem_type = mesh.hex_elements[elem_id].element_type
                hex_by_type.setdefault(elem_type, []).append(elem_id)

            for elem_type, elem_ids in hex_by_type.items():
                febio_type = FEBIO_HEX_TYPES.get(elem_type, "hex8")
                section_name = f"RingHex_{febio_type}"
                handler.write(f"    <Elements type=\"{febio_type}\" name=\"{section_name}\">\n")
                for element_id in sorted(elem_ids):
                    element = mesh.hex_elements[element_id]
                    conn = ",".join(str(nid) for nid in element.connectivity)
                    handler.write(f"      <elem id=\"{element.id}\">{conn}</elem>\n")
                handler.write("    </Elements>\n")

        if mesh.prism_elements:
            # Group by type
            prism_by_type: Dict[str, List[int]] = {}
            for elem_id in mesh.prism_elements:
                elem_type = mesh.prism_elements[elem_id].element_type
                prism_by_type.setdefault(elem_type, []).append(elem_id)

            for elem_type, elem_ids in prism_by_type.items():
                febio_type = FEBIO_PRISM_TYPES.get(elem_type, "penta6")
                section_name = f"RingPrism_{febio_type}"
                handler.write(f"    <Elements type=\"{febio_type}\" name=\"{section_name}\">\n")
                for element_id in sorted(elem_ids):
                    element = mesh.prism_elements[element_id]
                    conn = ",".join(str(nid) for nid in element.connectivity)
                    handler.write(f"      <elem id=\"{element.id}\">{conn}</elem>\n")
                handler.write("    </Elements>\n")

        # Surfaces
        if surfaces:
            for surf_name, facets in surfaces.items():
                handler.write(f"    <Surface name=\"{surf_name}\">\n")
                for elem_id, face_idx in facets:
                    # FEBio 4.0 surface facet format
                    handler.write(f"      <quad4 id=\"{elem_id}\">{elem_id}</quad4>\n")
                handler.write(f"    </Surface>\n")

        handler.write("  </Mesh>\n")
        
        # MeshDomains
        handler.write("  <MeshDomains>\n")
        if mesh.hex_elements:
            # Get unique element types
            hex_types = set(elem.element_type for elem in mesh.hex_elements.values())
            for elem_type in hex_types:
                febio_type = FEBIO_HEX_TYPES.get(elem_type, "hex8")
                section_name = f"RingHex_{febio_type}"
                handler.write(f"    <SolidDomain name=\"{section_name}\" mat=\"RingMaterial\"/>\n")
        if mesh.prism_elements:
            # Get unique element types
            prism_types = set(elem.element_type for elem in mesh.prism_elements.values())
            for elem_type in prism_types:
                febio_type = FEBIO_PRISM_TYPES.get(elem_type, "penta6")
                section_name = f"RingPrism_{febio_type}"
                handler.write(f"    <SolidDomain name=\"{section_name}\" mat=\"RingMaterial\"/>\n")
        handler.write("  </MeshDomains>\n")
        
        handler.write("</febio_spec>\n")
