# Inputs and software

## Input layout

Create the following directories in the repository root, or beneath `ADDOMER_WORKDIR` when using a separate work area:

```text
data/
├── 55aa_predictions.tar.gz
├── reference52_original/
│   ├── inputs.json
│   └── ... original files ending in seed_64438_sample_0.cif through sample_4.cif
└── reference_1X9T/
    ├── 1X9T.cif
    └── rcsb_pdb_1X9T.fasta
```

The archive is the original 55-aa download. Its expected members include `inputs.json`, `55aa/seed_79470/predictions/55aa_sample_0.cif` and `55aa_summary_confidence_sample_0.json`, repeated for samples 0–4. Step 01 writes normalised filenames to `results/predictions/` without altering the archive.

The laboratory reference is the **original** 52-aa ensemble, not a later reproduction. Preserve its original filenames so sample indices can be identified. Raw laboratory files are not bundled; possession of this code alone is insufficient to rerun the workflow.

Obtain the experimental structure from [RCSB 1X9T](https://www.rcsb.org/structure/1X9T). Use the original entry CIF containing reference chains A and B, not the large expanded biological assembly. The companion FASTA must contain a header and the extracted 523-residue reference sequence on its second line. This sequence starts `TGGRNS` and ends `SPRVLSSRTF`, corresponding to author residue numbers 49–571. It is not the complete deposited protein sequence. Step 02 asserts exact agreement with the saved alignment.

## Software environments

| Steps | Software |
|---|---|
| 01–06, 10–14 | Python 3.10 or newer and NumPy |
| 07 | Python environment providing PyMOL, `pymol2`, and `pymol.editor` |
| 08–09 | Python with NumPy and OpenMM 8.6.1; CPU platform |

A basic analysis environment can be installed with `python -m pip install -r requirements.txt`. OpenMM can be installed separately with `python -m pip install openmm==8.6.1`. Use an appropriate PyMOL distribution for step 07 and confirm `python -c "import pymol2; from pymol import editor"` succeeds in that environment. Installing NumPy does not install PyMOL.

OpenMM uses its bundled Amber14 and GBn2 parameters. Four CPU threads, seed 5500 for template 0, and bounded minimisation settings are retained. Software/platform differences may alter minimised coordinates. No model weights, MSA databases, Python installation or licensed software are redistributed.

## Outputs and reruns

All output is written below `results/`, with separate `predictions`, `interface`, `reconstruction`, `candidate_scan` and `candidate_geometry` directories. Existing outputs with the same names can be overwritten. Use a new work directory for a fresh analysis:

```bash
export ADDOMER_WORKDIR=/path/to/new_run
mkdir -p "$ADDOMER_WORKDIR/data"
# Place the input files there before executing scripts.
```

The Protenix server run itself is external to this repository. Its complete task record must be retained separately, including the six contact constraints that were absent from the downloaded input export. Input files, logs and raw structures should not be included in the code-only GitHub upload unless separately authorised for sharing.
