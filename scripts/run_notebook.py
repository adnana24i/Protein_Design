"""Execute notebooks/model_overview.ipynb in place so it ships with its outputs filled in.

Usage: python3 -I scripts/run_notebook.py
"""
import os
import sys

import nbformat
from nbclient import NotebookClient

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "notebooks", "model_overview.ipynb")

nb = nbformat.read(path, as_version=4)
NotebookClient(nb, timeout=900, kernel_name="python3",
               resources={"metadata": {"path": os.path.dirname(path)}}).execute()
nbformat.write(nb, path)

n_err = sum(1 for c in nb.cells for o in c.get("outputs", []) if o.get("output_type") == "error")
print(f"executed {path}: {len(nb.cells)} cells, {n_err} errors")
sys.exit(1 if n_err else 0)
