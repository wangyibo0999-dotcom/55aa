# Ordered execution guide

Run from the repository root. Use the Python environment specified for each stage. Scripts 01–06 and 10–14 use Python plus NumPy; 07 uses PyMOL Python; 08–09 use OpenMM Python. Steps operate on the original research inputs rather than arbitrary protein structures.

## 1 Audit prediction inputs

```bash
python scripts/01_audit_predictions.py
```

**Purpose:** verify sequence identities, summarise confidence and measure contacts in the five direct predictions. **Input:** the original 55-aa archive. **Output:** `results/predictions/sample_0.cif` through `sample_4.cif`, confidence JSONs, `inputs_exported.json` and `audit.json`. These normalised files feed the spatial checks.

## 2 Establish residue correspondence

```bash
python scripts/02_map_reference.py
```

**Purpose:** map 1X9T residues and test all penton chain permutations. **Input:** step 01 outputs, reference CIF/FASTA and original laboratory sample 0. **Output:** `results/interface/spatial_audit.json`. The original EMBOSS alignment is embedded and checked against the input sequences.

## 3 Check superposition and local sites

```bash
python scripts/03_check_superposition.py
```

**Purpose:** verify the rigid-fit implementation and compare local reference sites in all target and laboratory samples. **Input:** step 02 and all five original laboratory CIFs. **Output:** `results/interface/controls_and_local_fit.json`.

## 4 Scan the actual interfaces

```bash
python scripts/04_scan_interfaces.py
```

**Purpose:** enumerate actual contacts independently of the proposed reference contacts. **Input:** step 01 structures. **Output:** `results/interface/actual_interface.json`, including 4/6 Å contacts, contact consistency and atom pairs below 2 Å.

## 5 Generate reference position guides

```bash
python scripts/05_build_reference_guides.py
```

**Purpose:** transfer the reference peptide Cα positions onto a target penton. **Input:** steps 02–04. **Output:** `results/reconstruction/guide_manifest.json` and a penton PDB with five Cα-only guide chains. Those guide chains are not complete adapters.

## 6 Run the 300 rigid placements

```bash
python scripts/06_test_rigid_placements.py
```

**Purpose:** assess intact predicted trimers at the reference sites. **Input:** step 01 structures and step 05 guides. **Output:** `results/reconstruction/rigid_start_audit.json` and five diagnostic `guided55_source*_UNRELAXED.pdb` files.

Each of five trimers is fitted to 60 ordered assignments of three out of five sites: 300 location fits. Only the best three assignments per trimer receive detailed collision checks: 15 checks. These are rigid fits of existing coordinates, not 300 new predictions or flexible rebuilds. This branch did not produce accepted binding models.

## 7 Build complete target sequences from templates

```bash
# Run with Python configured for PyMOL and pymol2.
python scripts/07_build_template_models.py
```

**Purpose:** construct a separate template-derived hypothesis. **Input:** the five original laboratory 52-aa structures. **Output:** five `results/reconstruction/template55_source*_UNRELAXED.pdb` files and `template_build.json`. K27L, N45D and T52A substitutions plus terminal DNA extension are applied to every adapter; sequences are checked. This stage does not use the rigid placements from step 06.

## 8 Perform the first local minimisation

```bash
# Run with OpenMM 8.6.1 and NumPy.
python scripts/08_minimise_fixed_core.py 0
```

**Purpose:** reduce local strain with much of the scaffold fixed. **Input:** unrelaxed template 0. **Output:** `template55_source0_prepared.pdb`, `template55_source0_LOCAL_MINIMIZED.pdb` and `relaxation_0.json` under `results/reconstruction/`. The minimisation is bounded at 200 iterations.

## 9 Allow penton amino termini to move

```bash
python scripts/09_minimise_flexible_termini.py 0
```

**Purpose:** address inherited terminal overlaps with penton residues 1–48 mobile. **Input:** step 08 `LOCAL_MINIMIZED.pdb`, consumed directly. **Output:** `template55_source0_FLEX_NTERM_MINIMIZED.pdb` and `relaxation_flex_0.json`. The minimisation is bounded at 500 iterations. Completion is not proof of convergence or binding validity.

## 10 Check reconstruction integrity

```bash
python scripts/10_check_rebuilt_models.py
```

**Purpose:** check sequences, peptide connections, reference contact distances and adapter–penton short contacts. **Input:** reconstructed PDBs and step 01 input sequences. **Output:** `results/reconstruction/rebuilt_model_checks.json`.

## 11 Audit all interchain overlaps

```bash
python scripts/11_audit_interchain_clashes.py
```

**Purpose:** include penton–penton and adapter–adapter overlaps as well as adapter–penton overlaps. **Input:** prepared, first-stage and final template 0 coordinates. **Output:** `results/reconstruction/all_interchain_short_contacts.json`.

## 12 Enumerate candidate residue pairs

```bash
python scripts/12_scan_candidate_pairs.py
```

**Purpose:** identify hypotheses from contacts in the original laboratory templates and check sequence correspondence to the target. **Input:** original laboratory CIFs and `inputs.json`. **Output:** `results/candidate_scan/scan.json` and `rank.json`. This preliminary screen uses a 10° grid and 2 Å sulfur clearance; it is not the final target check.

## 13 Refine preliminary pairs

```bash
python scripts/13_refine_candidate_pairs.py
```

**Purpose:** re-evaluate initially supported pairs with a 5° grid and 2.5 Å sulfur clearance. **Input:** step 12 outputs and original laboratory CIFs. **Output:** `results/candidate_scan/refined_stricter.json`. Its rankings are geometric observations, not probabilities of stability improvement.

## 14 Check the two target candidates

```bash
python scripts/14_check_final_candidates.py
```

**Purpose:** assess V31C–Q457C and E30C–L223C in all 15 adapter–penton chain combinations of each target model. **Input:** five unrelaxed templates from 07 and the final minimised representative from 09. **Output:** `results/candidate_geometry/audit.json` and passing single-chain-pair mutant PDBs.

The final screen uses a 2° grid, sulfur separation 1.8–2.3 Å, both bond angles 85–125°, χ3 within 30° of ±90°, and at least 2.5 Å sulfur clearance. Passing PDBs are inspection structures, without mutant-specific relaxation or evidence of successful covalent linkage. See the export correction in `METHODS_AND_LIMITATIONS.md` before comparing these PDBs with older inspection files.
