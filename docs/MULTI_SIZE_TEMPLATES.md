# Multi-Size 2D Mesh Templates — Implementation Guide

This document describes how to implement generation of 2D meshed templates for **all implant sizes** from a parameter table, and how to export them for use in another program. Use it when adding multi-size and batch-export support to the templateator pipeline.

---

## 1. Context

- **CAD**: Fusion 360 designs are fully parametrized; sizes are controlled via a **parameter table** (e.g. User Parameters).
- **Fixed sizes**: There is a finite, known set of sizes (e.g. 148, 150, 152, …).
- **Goal**: Produce one `.json` template per size (`template_{id}.json`) and make them available to another program (e.g. by copying into that program’s repo).

**Current pipeline (single size):**

1. Fusion 360: intersect body with `midplane_real` → export `boundary.json`
2. `mesh_with_gmsh.py`: `boundary.json` → `template.json` (nodes, elements, surfaces)
3. Downstream: reads `template.json`

**Target pipeline (all sizes):**

1. **Parameter table** defines each size and, for Fusion, the parameter overrides.
2. **Boundaries**: one `ring_{id}_boundary.json` per size (from Fusion, or from a parametric generator).
3. **Batch mesh**: run GMSH for each boundary → `template_{id}.json` + `index.json`.
4. **Export**: copy `template_*.json` and `index.json` into the other program’s repo.

---

## 2. Parameter Table (`config/sizes.json`)

Central list of sizes and parameters. Format is meant to support both **Fusion-driven** and (optionally) **parametric boundary generation**.

### 2.1 Structure

```json
{
  "description": "Parameter table for medical implant sizes. Adapt fusion_params keys to match Fusion 360 User Parameter names.",
  "sizes": [
    {
      "id": "148",
      "fusion_params": {
        "Size": 148,
        "Thickness": 1.9
      }
    },
    {
      "id": "150",
      "fusion_params": {
        "Size": 150,
        "Thickness": 2.0
      }
    }
  ],
  "mesh": {
    "divisions": 5,
    "element_size": null,
    "quality": false
  }
}
```

- **`id`**: Unique size identifier; used in filenames: `ring_{id}_boundary.json`, `template_{id}.json`.
- **`fusion_params`**: For Fusion 360: `UserParameter` **names** → values. These drive the CAD when doing batch export. Keys must match the names in Fusion’s Parameters panel.
- **`mesh`**: Defaults for `mesh_with_gmsh` (divisions, element_size, quality). Scripts can override via CLI.

Optional extension for a **parametric boundary generator** (no Fusion):

- **`geometry_params`**: e.g. `{ "inner_r": 7.5, "outer_r": 9.0, "thickness": 2.0, "fillet_r": 0.5 }` if the 2D section can be fully defined by formulas. Used only when generating `boundary.json` from geometry, not by Fusion.

---

## 3. File Naming and Layout

| Role            | Path / pattern                         | Produced by                      |
|-----------------|----------------------------------------|----------------------------------|
| Parameter table | `config/sizes.json`                    | Hand-authored / from design DB   |
| Boundary (Fusion)| `boundaries/ring_{id}_boundary.json`  | Fusion 360 export (single/batch) |
| Legacy boundary | `boundaries/boundary.json`             | Fusion single export (old name)  |
| Mesh template   | `output_templates/template_{id}.json`  | `mesh_with_gmsh` (batch script)  |
| Index           | `output_templates/index.json`          | Batch mesh script                |

- **Index** (`index.json`): `{ "sizes": ["148","150",…], "templates": { "148": "template_148.json", … }, "source": "templateator batch_mesh_all_sizes" }` so the other program can discover which sizes exist and which file to load.

---

## 4. Fusion 360: Boundary Export (Single and Batch)

### 4.1 Single Export (Current Behavior)

- User sets size in Fusion (e.g. by changing parameters or loading a configuration).
- Run the add-in → intersect body with `midplane_real` → save boundary via “Save” dialog.
- User can save as `ring_{id}_boundary.json` in `boundaries/` to align with multi-size naming.

No code changes required for single export; only the chosen filename should follow `ring_{id}_boundary.json` when preparing for batch.

### 4.2 Batch Export (To Be Implemented)

Idea: one Fusion run that exports **all** sizes in `config/sizes.json` without reopening the document.

**Prerequisites:**

- `SCRIPT_DIR` and `BOUNDARIES_FOLDER`:
  - `SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))`
  - `BOUNDARIES_FOLDER = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "boundaries"))`
- Design is parametric; `fusion_params` keys match **User Parameter** names (as in the Parameters panel).
- Construction plane `midplane_real` and at least one body exist.

**Algorithm (run_batch):**

1. Load `config/sizes.json` (e.g. `config_path = os.path.join(SCRIPT_DIR, "..", "config", "sizes.json")`).
2. `os.makedirs(BOUNDARIES_FOLDER, exist_ok=True)`.
3. For each `s` in `sizes` with `s["id"]`:
   - **Set parameters**: for each `(k, v)` in `s.get("fusion_params", {})`:  
     `p = design.userParameters.itemByName(k)`; if `p`: `p.expression = str(v)`.
   - **Regenerate**: run Fusion’s Regenerate (e.g. `commandDefinitions.itemById("FusionRegenerateCmd")` or `"DesignRegenerate"` and `execute()`). If the API does not expose it reliably, document that the user may need to ensure the design is up to date (e.g. one manual Regen before starting).
   - **Body**: use e.g. `root_comp.bRepBodies.item(0)`. Optionally support `body_index` or `body_name` in config later.
   - **Extract boundary**: same steps as single export (intersect with `midplane_real`, outer profile, `extract_boundary_curves` / `extract_boundary_vertices`, `build_boundary_json`). The only difference: pass `name=f"ring_{id}_boundary"` into `build_boundary_json` so the `"name"` field in JSON is correct.
   - **Write**: `boundaries/ring_{id}_boundary.json` (no Save dialog). Delete the temporary intersection sketch after each size.
4. Show a summary: exported count and list of errors (e.g. missing parameter, no profile, no body).

**Refactor to support batch:**

- `build_boundary_json(curves, vertices, name=None)`: if `name` is `None`, use `TEMPLATE_NAME`; otherwise use `name`. Use it for both the top-level `"name"` and in `metadata.description`.
- `export_boundary_data(design, template_name, body=None, body_index=0)`:
  - Find `midplane_real`, resolve `body` (or `bodies.item(body_index)`), create intersection sketch, get outer profile, extract curves/vertices, call `build_boundary_json(..., name=template_name)`.
  - Return `(data_dict, error_str, sketch, created_new)`. On failure: `(None, message, sketch_or_None, created)`. Do **not** write to disk; do **not** delete the sketch. Caller deletes the sketch and writes the file.
- `run()` (single): after resolving `target_body`, call `export_boundary_data(design, TEMPLATE_NAME, body=target_body)`. On success: `get_output_path(ui)`, write JSON, summary message, then “Delete temporary sketch?” and optional `sketch.deleteMe()`.
- `run_batch(design, ui, config_path)`: implement as above; use `export_boundary_data(design, f"ring_{id}_boundary", body=body)` and write to `BOUNDARIES_FOLDER`.
- **Batch entry point**: if `config/sizes.json` exists and has more than one size with an `id`, show a dialog: “Run batch export for N sizes? (No = single export for current design)”. If Yes → `run_batch(design, ui, config_path)` and return; else continue with the existing single-export flow.
- **Save dialog / default directory**: in `get_output_path(ui)`, when not in batch, set `dialog.initialDirectory = BOUNDARIES_FOLDER` (and `os.makedirs(BOUNDARIES_FOLDER, exist_ok=True)` if needed) so the default save location is `boundaries/`.

**Fusion API notes:**

- `design.userParameters.itemByName("Size")` and `p.expression = "150"` (string). `itemById` is an alternative if you use internal IDs.
- Regenerate: command IDs can differ by Fusion version; a try/except around `execute()` is recommended.

---

## 5. Batch Mesh Script (`scripts/batch_mesh_all_sizes.py`)

**Role:** For each size in `config/sizes.json`, find `boundaries/ring_{id}_boundary.json`, run `mesh_with_gmsh.generate_mesh(...)`, write `output_templates/template_{id}.json`. At the end, write `output_templates/index.json`.

**Inputs:**

- `--config`: path to `sizes.json` (default: `config/sizes.json`).
- `--boundaries-dir`: where to find `ring_{id}_boundary.json` (default: `boundaries/`).
- `--output-dir`: where to write `template_{id}.json` and `index.json` (default: `output_templates/`).
- `-n` / `--divisions`, `--element-size`, `-q` / `--quality`, `--no-quality`: override `config.mesh` for GMSH.

**Logic:**

- Load config, read `sizes` and `mesh` defaults.
- For each `s` with `s["id"]`:
  - `boundary_file = boundaries_dir / f"ring_{sid}_boundary.json"`.
  - Optional fallback: if that file is missing, `len(sizes)==1`, and `boundaries_dir / "boundary.json"` exists, use `boundary.json` (backward compatibility for a single-size setup with the legacy name).
  - If the chosen boundary file still does not exist: **skip** that size, print `[skip] {id}: boundary not found`.
  - Else: `generate_mesh(boundary_file, output_dir / f"template_{id}.json", divisions=..., element_size=..., quality=..., preview=False)`. On exception: collect error, **skip**.
- Write `output_dir / "index.json"` with `sizes`, `templates` map, and a `source` string.
- Print summary: generated list, skipped list, errors. Exit with non-zero if there are errors.

**Dependency:** `from mesh_with_gmsh import generate_mesh`. Ensure `scripts/` is on `sys.path` when running the script (e.g. `sys.path.insert(0, str(SCRIPT_DIR))` at the top).

---

## 6. Export to Other Repo (`scripts/export_to_repo.py`)

**Role:** Copy `template_*.json` and optionally `index.json` from `output_templates/` to a target directory (e.g. the other program’s `templates/` or `assets/`).

**Usage:**

```bash
python scripts/export_to_repo.py --target /path/to/other/repo/templates
python scripts/export_to_repo.py --target ../downstream_app/templates --no-index
```

**Arguments:**

- `--target` / `-t`: required; destination directory. Create with `mkdir(parents=True, exist_ok=True)`.
- `--source` / `-s`: source folder (default: `output_templates/`).
- `--include-index` (default: True): copy `index.json`.
- `--no-index`: do not copy `index.json`.
- `--glob`: pattern for template files (default: `template_*.json`).

**Behavior:** Copy each matching file and (if enabled) `index.json` with `shutil.copy2`. Print list of copied files. No parsing of JSON.

---

## 7. End-to-End Workflow

**Option A — Fusion batch + batch mesh + export**

1. In Fusion: open the parametrized design, run the add-in. If `config/sizes.json` has multiple sizes, choose “Yes” for batch → all `ring_{id}_boundary.json` are written to `boundaries/`.
2. In a terminal:  
   `python scripts/batch_mesh_all_sizes.py`  
   → `output_templates/template_{id}.json` and `index.json`.
3. Export:  
   `python scripts/export_to_repo.py --target /path/to/other/repo/templates`

**Option B — Single-size Fusion + batch mesh**

1. For each size: in Fusion, set parameters, run add-in, save as `boundaries/ring_{id}_boundary.json`.
2. Run `batch_mesh_all_sizes.py` and then `export_to_repo.py` as above.

**Option C — Parametric boundaries (no Fusion) — future**

- If the 2D section is expressible as a function of `geometry_params` (e.g. arcs, lines, fillets), add `scripts/generate_boundary_from_params.py` that, for each size, generates `boundaries/ring_{id}_boundary.json` in the same JSON schema as the Fusion export. Then the rest of the pipeline (batch mesh, export) is unchanged.

---

## 8. What the Other Program Needs

- **Files:** `template_{id}.json` (and optionally `index.json`).
- **Template format:** existing `template.json` schema: `name`, `version`, `units`, `profile_plane`, `base_height_mm`, `metadata`, `nodes`, `elements`, `boundary_nodes`, `surfaces`. `ring_template_reader_json.py` (or equivalent) can stay as is; it only needs the path to the right `template_{id}.json`.
- **Discovery:** if `index.json` is exported, the other program can read `sizes` and `templates` to know which sizes exist and which filename to use for each.

---

## 9. Summary of New / Modified Artifacts

| Artifact               | Action                                      |
|------------------------|---------------------------------------------|
| `config/sizes.json`    | **Create**: parameter table + mesh defaults |
| `scripts/batch_mesh_all_sizes.py` | **Create**: batch GMSH + index       |
| `scripts/export_to_repo.py`       | **Create**: copy to target directory |
| `fusion360_ring_template_export/fusion360_ring_template_export.py` | **Modify**: `build_boundary_json(name=...)`, `export_boundary_data(...)`, `run_batch(...)`, batch prompt, `get_output_path` default dir, `SCRIPT_DIR`/`BOUNDARIES_FOLDER` |

---

## 10. Optional: Parametric Boundary Generator

If the cross-section can be fully defined by a few numbers (e.g. inner/outer radius, thickness, fillet radius) that are in the parameter table, a **headless** path is possible:

- Add `geometry_params` to each size in `sizes.json` (or derive from `fusion_params` via a separate lookup table).
- **`scripts/generate_boundary_from_params.py`**: for each size, build a `boundary` dict (same structure as Fusion: `curves` with `line`/`arc`/`spline`, `vertices`) from the geometry params, and write `boundaries/ring_{id}_boundary.json`. The exact formulas depend on the implant’s 2D profile (e.g. trapezoid with rounded corners, custom splines).
- Then run `batch_mesh_all_sizes.py` and `export_to_repo.py` as usual. This avoids opening Fusion for each size and can be run in CI.

This is only feasible when the shape is formula-driven; it is not described in full here.
