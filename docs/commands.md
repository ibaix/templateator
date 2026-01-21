## Workflow

1. **Fusion 360**: Run the export script to intersect body with midplane_real → exports boundary.json

2. **Generate mesh**:
```bash
.\venv\Scripts\python.exe scripts/mesh_with_gmsh.py boundaries/boundary.json -o output_templates/template.json
```

3. **Visualize**:
```bash
.\venv\Scripts\python.exe scripts/visualize_template.py output_templates/template_2.json
```

## Mesh Generation Options

```bash
# Basic mesh with 5 divisions
python scripts/mesh_with_gmsh.py boundary.json -o template.json -n 5

# Preview mesh in GMSH GUI
python scripts/mesh_with_gmsh.py boundary.json -o template.json --preview

# HIGH QUALITY mesh (recommended for FEM simulations)
python scripts/mesh_with_gmsh.py boundary.json -o template.json -q --min-quality 0.5

# Fine mesh with quality optimization
python scripts/mesh_with_gmsh.py boundary.json -o template.json -n 10 -q --min-quality 0.4 --optimize-iter 5
```

## Quality Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `-q, --quality` | off | Enable Laplace + Relocate smoothing |
| `--min-quality` | 0.3 | Min quality for quad recombination (0-1). Higher = fewer degenerate quads |
| `--optimize-iter` | 3 | Number of smoothing iterations |

## Quality Thresholds

For FEM simulations, aim for:
- **Angles**: 60°-120° (ideal 90°)
- **Quality score**: > 0.5
- **No elements with angles >150° or <30°**