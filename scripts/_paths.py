"""Project-relative paths. Optional ADDOMER_WORKDIR relocates data/results only."""
from pathlib import Path
import os
SCRIPTS = Path(__file__).resolve().parent
REPOSITORY = SCRIPTS.parent
WORK = Path(os.environ.get("ADDOMER_WORKDIR", REPOSITORY)).resolve()
DATA = WORK / "data"
RESULTS = WORK / "results"
def output_dir(name):
    path = RESULTS / name
    path.mkdir(parents=True, exist_ok=True)
    return path
