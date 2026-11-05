"""Re-run the finite joint sector search and the dense broadband audit."""
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

from root_experiments import coverage

if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
    out = root / "outputs" / "coverage" / stamp
    out.mkdir(parents=True)
    coverage(out)
    print(out)
