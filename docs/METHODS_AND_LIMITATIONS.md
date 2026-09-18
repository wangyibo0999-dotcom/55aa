# Methods and interpretation limits

## Core scope

This export covers the revised 55-aa analysis: input checks, reference mapping, interface assessment, rigid-placement diagnostics, template reconstruction, bounded minimisation and candidate screening. The original 52-aa ensemble remains necessary because it supplied comparison structures and reconstruction templates.

The following are intentionally omitted: historical 50-aa PyRosetta/Slurm experiments, duplicated scheduler configurations, repeated 52-aa server-run comparisons, intermediate reporting scripts, PyMOL presentation scripts, manuscript generation and private raw data. Those records remain in the original project archive. This compact repository is not a complete archive of every calculation mentioned in the report.

## Scientific interpretation

- 1X9T provides homologous peptide-position evidence, not an experimentally determined full 55-aa complex.
- The direct predictions and reconstructed coordinates are separate evidence sources.
- The 300 placements comprise rigid fitting trials; only 15 received detailed collision checks. The diagnostic placements were not accepted as binding models.
- Reconstruction modifies inherited template coordinates. Only template 0 underwent the two bounded minimisation stages. The calculation is not molecular dynamics or binding free-energy estimation.
- Atom pairs below 2 Å flag severe short contacts; their absence is not full stereochemical validation.
- Candidate sulfur geometry is an operational screen, not an experimental success probability or stability estimate. The minimised representative and candidate-supporting templates are different conformations.
- Multiple copies in a template are not independent experimental replicates. No numerical probability of experimental success was established.

## Changes made during code curation

1. Files were renamed by execution order. Personal absolute paths were replaced with `data/` and `results/` paths managed by `_paths.py`. The scientific source file and original hash are recorded in `metadata/CODE_PROVENANCE.json`.
2. The first OpenMM output is consumed directly by step 09. The earlier workflow manually copied the same bytes from a `LOCAL_MINIMIZED` filename to `FIXED_CORE_TRIAL`; that naming-only handoff is no longer needed.
3. The final candidate check consumes `FLEX_NTERM_MINIMIZED.pdb` directly rather than a separately copied delivery file. Its numerical identity to the original delivery structure was checked during packaging.
4. Only the atom and geometry helper definitions required by these scripts were retained from the larger historical screen. Narrative delivery/report generation was omitted where it was not needed for downstream calculations.
5. Export starts from the original working scripts. A previous export's broad email substitution had corrupted an `@` matrix multiplication expression. No blanket text redaction is applied to operators in this version, and step 06 was executed on the original inputs.

## Final geometry threshold correction

The earlier final-check script attempted to change the sulfur clearance from 2.0 to 2.5 Å by replacing `md>=2` in the output of `ast.unparse`. That output contains spaces (`md >= 2`), so the replacement did not take effect. Step 14 now matches the actual expression and asserts that every passing pose meets 2.5 Å clearance.

The corrected script was rerun against the same six saved structures. Candidate pass counts were unchanged: V31C–Q457C passed three pairs in template 2, and E30C–L223C passed three in template 4; both had zero passes in the minimised representative. One previously exported E30C–L223C pose in template 4, chains G–C, had clearance approximately 2.475 Å. The corrected search selects a different qualifying pose. Consequently, the pass-count conclusion is unchanged, but old inspection coordinates must not be treated as the corrected outputs. Generate current inspection PDBs with step 14.

## Validation scope

The export was syntax-checked and steps 01–06 and 10–14 were executed in an isolated work area using the original study inputs and saved reconstruction coordinates. Steps 07–09 were not rerun during packaging; their original calculations are distinct from export validation. The full Protenix server prediction, PyMOL reconstruction and OpenMM minimisation have therefore not been repeated end-to-end for this export. See `metadata/VALIDATION.json` for recorded checks.
