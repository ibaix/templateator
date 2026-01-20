"""
Ring Template Visualizer

Visualizes the mesh template showing nodes and elements (quads and triangles) in a 2D plot.

Usage (from project root):
    python scripts/visualize_template.py output_templates/template.json
    python scripts/visualize_template.py output_templates/template.json --output visualizations/

Options:
    --output, -o    Output directory for visualization PNG (default: same as input)
    --no-ids        Hide node ID labels
    --show-elem-ids Show element ID labels
"""

import json
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np


def load_template(filepath):
    """Load template from JSON file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        return json.load(f)


def visualize_template(template, show_node_ids=True, show_element_ids=False):
    """
    Visualize the mesh template.
    
    Parameters
    ----------
    template : dict
        The loaded JSON template
    show_node_ids : bool
        Whether to display node ID labels
    show_element_ids : bool
        Whether to display element ID labels at centroids
    """
    nodes = template.get('nodes', [])
    elements = template.get('elements', [])
    
    # Create node lookup dictionary
    node_map = {n['id']: (n['x'], n['z']) for n in nodes}
    
    # Set up the figure
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    
    # Collect patches for quads and tris
    quad_patches = []
    tri_patches = []
    
    quad_count = 0
    tri_count = 0
    
    for elem in elements:
        elem_type = elem['type']
        node_ids = elem['nodes']
        
        # Get coordinates for each node in the element
        try:
            coords = [node_map[nid] for nid in node_ids]
        except KeyError as e:
            print(f"Warning: Element {elem['id']} references missing node {e}")
            continue
        
        # Create polygon
        polygon = patches.Polygon(coords, closed=True)
        
        if elem_type == 'QUAD':
            quad_patches.append(polygon)
            quad_count += 1
        elif elem_type == 'TRI':
            tri_patches.append(polygon)
            tri_count += 1
        
        # Show element ID at centroid
        if show_element_ids:
            centroid_x = np.mean([c[0] for c in coords])
            centroid_z = np.mean([c[1] for c in coords])
            ax.text(centroid_x, centroid_z, str(elem['id']), 
                   fontsize=6, ha='center', va='center', color='black', alpha=0.7)
    
    # Add quad patches with blue fill
    if quad_patches:
        quad_collection = PatchCollection(quad_patches, 
                                          facecolor='#4A90D9', 
                                          edgecolor='#1a1a2e', 
                                          linewidth=0.8,
                                          alpha=0.6)
        ax.add_collection(quad_collection)
    
    # Add tri patches with orange fill
    if tri_patches:
        tri_collection = PatchCollection(tri_patches, 
                                         facecolor='#E87D3E', 
                                         edgecolor='#1a1a2e', 
                                         linewidth=0.8,
                                         alpha=0.6)
        ax.add_collection(tri_collection)
    
    # Plot nodes
    node_x = [n['x'] for n in nodes]
    node_z = [n['z'] for n in nodes]
    ax.scatter(node_x, node_z, c='#2d2d44', s=15, zorder=5, edgecolor='white', linewidth=0.5)
    
    # Show node IDs
    if show_node_ids:
        for n in nodes:
            ax.annotate(str(n['id']), 
                       (n['x'], n['z']), 
                       fontsize=5, 
                       ha='left', 
                       va='bottom',
                       xytext=(1, 1),
                       textcoords='offset points',
                       color='#1a1a2e',
                       alpha=0.8)
    
    # Set equal aspect ratio and labels
    ax.set_aspect('equal')
    ax.set_xlabel('X (mm)', fontsize=12)
    ax.set_ylabel('Z (mm)', fontsize=12)
    
    # Title with statistics
    title = f"{template.get('name', 'Template')}\n"
    title += f"Nodes: {len(nodes)} | Elements: {len(elements)} (Quads: {quad_count}, Tris: {tri_count})"
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        patches.Patch(facecolor='#4A90D9', edgecolor='#1a1a2e', alpha=0.6, label=f'QUAD ({quad_count})'),
        patches.Patch(facecolor='#E87D3E', edgecolor='#1a1a2e', alpha=0.6, label=f'TRI ({tri_count})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2d2d44', 
               markersize=8, label=f'Nodes ({len(nodes)})')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Grid
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig, ax


def main():
    import argparse
    from pathlib import Path
    
    parser = argparse.ArgumentParser(
        description="Visualize mesh template showing nodes and elements",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (from project root):
    python scripts/visualize_template.py output_templates/template.json
    python scripts/visualize_template.py output_templates/template.json -o visualizations/
    python scripts/visualize_template.py output_templates/template.json --no-ids
        """
    )
    parser.add_argument("template_file", nargs='?', default='output_templates/template.json',
                        help="Input template JSON file (default: output_templates/template.json)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output directory for PNG (default: same directory as input)")
    parser.add_argument("--no-ids", action="store_true",
                        help="Hide node ID labels")
    parser.add_argument("--show-elem-ids", action="store_true",
                        help="Show element ID labels")
    parser.add_argument("--no-show", action="store_true",
                        help="Don't show interactive plot (just save PNG)")
    
    args = parser.parse_args()
    template_file = args.template_file
    
    print(f"Loading template: {template_file}")
    
    try:
        template = load_template(template_file)
    except FileNotFoundError:
        print(f"Error: File '{template_file}' not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in '{template_file}': {e}")
        sys.exit(1)
    
    # Print summary
    nodes = template.get('nodes', [])
    elements = template.get('elements', [])
    quad_count = sum(1 for e in elements if e['type'] == 'QUAD')
    tri_count = sum(1 for e in elements if e['type'] == 'TRI')
    
    print(f"\nTemplate Summary:")
    print(f"  Name: {template.get('name', 'N/A')}")
    print(f"  Nodes: {len(nodes)}")
    print(f"  Elements: {len(elements)}")
    print(f"    - QUADs: {quad_count}")
    print(f"    - TRIs: {tri_count}")
    
    # Visualize
    show_node_ids = not args.no_ids
    fig, ax = visualize_template(template, show_node_ids=show_node_ids, show_element_ids=args.show_elem_ids)
    
    # Determine output path
    input_path = Path(template_file)
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{input_path.stem}_visualization.png"
    else:
        output_file = input_path.with_suffix('.png').with_stem(input_path.stem + '_visualization')
    
    fig.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nVisualization saved to: {output_file}")
    
    # Show interactive plot
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
