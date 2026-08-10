import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import similarityengine as sim  # noqa: E402

print("=== Signature de find_similar_genes ===")
print(inspect.signature(sim.find_similar_genes))

print("\n=== Code source de find_similar_genes ===")
print(inspect.getsource(sim.find_similar_genes))

print("\n=== Signature de compare_with_database ===")
print(inspect.signature(sim.compare_with_database))
