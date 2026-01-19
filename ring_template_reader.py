"""Herramientas para leer plantillas *.inp de Abaqus."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Tuple


class InpFormatError(RuntimeError):
    """Excepción personalizada para errores de formato en archivos .inp."""


@dataclass(frozen=True)
class NodeRecord:
    """Representa un nodo leído desde la plantilla."""

    id: int
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class ElementRecord:
    """Representa un elemento shell proveniente de la plantilla."""

    id: int
    connectivity: Tuple[int, ...]
    elset: str | None = None


@dataclass(frozen=True)
class MeshTemplate:
    """Agrupa los nodos y elementos extraídos de la plantilla."""

    nodes: Dict[int, NodeRecord]
    elements: Dict[int, ElementRecord]
    element_type: str
    surfaces: Dict[str, Tuple[Tuple[int, int], ...]] = field(default_factory=dict)
    element_order: int = 1  # 1=linear (QUAD4/TRI3), 2=quadratic (QUAD8/TRI6)


def _parse_option(line: str, key: str) -> str | None:
    """Busca el valor de una opción en la línea de cabecera."""

    key_upper = key.upper()
    for option in line.split(","):
        option = option.strip()
        if "=" not in option:
            continue
        k, v = option.split("=", 1)
        if k.strip().upper() == key_upper:
            return v.strip()
    return None


def _ensure_integer(value: float) -> int:
    """Verifica que un valor flotante represente un entero exacto."""

    if int(value) != value:
        raise InpFormatError(
            f"Los identificadores deben ser enteros, se recibió {value!r}."
        )
    return int(value)


def _parse_node_line(line: str) -> NodeRecord:
    parts = [component.strip() for component in line.split(",")]
    if len(parts) < 4:
        raise InpFormatError(
            "Cada línea de nodos debe contener ID seguido de tres coordenadas."
        )
    node_id = _ensure_integer(float(parts[0]))
    x, y, z = (float(parts[i]) for i in range(1, 4))
    return NodeRecord(node_id, x, y, z)


def _parse_element_line(line: str, expected_nodes: int) -> Tuple[int, Tuple[int, ...]]:
    parts = [component.strip() for component in line.split(",") if component.strip()]
    if len(parts) < expected_nodes + 1:
        raise InpFormatError(
            "Cada elemento debe incluir un ID seguido de su conectividad completa."
        )
    element_id = _ensure_integer(float(parts[0]))
    try:
        connectivity = tuple(int(node) for node in parts[1:])
    except ValueError as exc:  # pragma: no cover - manejo defensivo
        raise InpFormatError("La conectividad del elemento debe ser numérica.") from exc

    if len(connectivity) != expected_nodes:
        raise InpFormatError(
            f"Se esperaba conectividad de {expected_nodes} nodos y se obtuvieron {len(connectivity)}."
        )

    return element_id, connectivity


def read_inp(path: str | Path, *, shell_nodes: int = 4) -> MeshTemplate:
    """Lee una plantilla `.inp` y devuelve los nodos y elementos shell.

    Parameters
    ----------
    path:
        Ruta al archivo de plantilla.
    shell_nodes:
        Número de nodos esperados para cada elemento shell.
    """

    filepath = Path(path)
    if not filepath.exists():
        raise FileNotFoundError(f"No se encontró la plantilla: {filepath}")

    nodes: Dict[int, NodeRecord] = {}
    elements: Dict[int, ElementRecord] = {}
    element_type: str | None = None
    element_order: int = 1  # 1=linear, 2=quadratic
    current_elset: str | None = None
    in_node_block = False
    in_element_block = False
    in_surface_block = False
    current_surface: str | None = None
    surfaces: Dict[str, list[Tuple[int, int]]] = {}
    current_shell_nodes = shell_nodes

    with filepath.open("r", encoding="utf-8") as handler:
        for raw_line in handler:
            stripped = raw_line.strip()

            if not stripped:
                continue
            if stripped.startswith("**"):
                continue

            upper = stripped.upper()
            if upper.startswith("*NODE"):
                in_node_block = True
                in_element_block = False
                in_surface_block = False
                continue
            if upper.startswith("*ELEMENT"):
                in_node_block = False
                in_element_block = True
                in_surface_block = False
                element_type = _parse_option(stripped, "TYPE") or element_type
                current_elset = _parse_option(stripped, "ELSET")
                # Detect shell node count and order from element type
                type_option = _parse_option(stripped, "TYPE")
                if type_option:
                    type_upper = type_option.upper()
                    # Second-order (quadratic) elements
                    if "S8" in type_upper or "CPS8" in type_upper or "CPE8" in type_upper:
                        current_shell_nodes = 8
                        element_order = 2
                    elif "S6" in type_upper or "CPS6" in type_upper or "CPE6" in type_upper:
                        current_shell_nodes = 6
                        element_order = 2
                    # First-order (linear) elements
                    elif "S4" in type_upper or "CPS4" in type_upper or "CPE4" in type_upper:
                        current_shell_nodes = 4
                        element_order = 1
                    elif "S3" in type_upper or "CPS3" in type_upper or "CPE3" in type_upper:
                        current_shell_nodes = 3
                        element_order = 1
                    else:
                        current_shell_nodes = shell_nodes
                continue
            if upper.startswith("*SURFACE"):
                in_node_block = False
                in_element_block = False
                in_surface_block = True
                current_surface = _parse_option(stripped, "NAME")
                if current_surface:
                    surfaces.setdefault(current_surface, [])
                continue
            if upper.startswith("*"):
                in_node_block = False
                in_element_block = False
                in_surface_block = False
                continue

            if in_node_block:
                node = _parse_node_line(stripped)
                if node.id in nodes:
                    raise InpFormatError(f"ID de nodo duplicado encontrado: {node.id}")
                nodes[node.id] = node
                continue

            if in_element_block:
                element_id, connectivity = _parse_element_line(stripped, current_shell_nodes)
                if element_id in elements:
                    raise InpFormatError(f"ID de elemento duplicado encontrado: {element_id}")
                elements[element_id] = ElementRecord(
                    element_id, connectivity, elset=current_elset
                )
                continue

            if in_surface_block and current_surface:
                parts = [p.strip() for p in stripped.split(",") if p.strip()]
                if len(parts) != 2 or not parts[1].upper().startswith("S"):
                    # Skip malformed surface lines
                    continue
                elem_id = _ensure_integer(float(parts[0]))
                try:
                    face_idx = int(parts[1][1:])
                except ValueError:
                    continue
                surfaces[current_surface].append((elem_id, face_idx))
                continue

    if not nodes:
        raise InpFormatError("El archivo no contiene un bloque *NODE válido.")
    if not elements:
        raise InpFormatError("El archivo no contiene un bloque *ELEMENT válido.")
    if element_type is None:
        raise InpFormatError("No se pudo determinar el tipo de elemento de la plantilla.")

    # Normalize surfaces to tuples
    normalized_surfaces: Dict[str, Tuple[Tuple[int, int], ...]] = {
        name: tuple(entries) for name, entries in surfaces.items()
    }

    return MeshTemplate(
        nodes=nodes,
        elements=elements,
        element_type=element_type,
        surfaces=normalized_surfaces,
        element_order=element_order,
    )


def max_node_id(nodes: Iterable[NodeRecord]) -> int:
    """Obtiene el máximo identificador de nodo de una colección."""

    try:
        return max(node.id for node in nodes)
    except ValueError as exc:  # pragma: no cover - colecciones vacías ya fueron validadas
        raise InpFormatError("La plantilla no contiene nodos.") from exc


def max_element_id(elements: Iterable[ElementRecord]) -> int:
    """Obtiene el máximo identificador de elemento de una colección."""

    try:
        return max(element.id for element in elements)
    except ValueError as exc:  # pragma: no cover - colecciones vacías ya fueron validadas
        raise InpFormatError("La plantilla no contiene elementos.") from exc

