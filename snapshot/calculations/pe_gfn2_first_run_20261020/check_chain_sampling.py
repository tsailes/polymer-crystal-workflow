"""A limited Gamma-supercell sampling diagnostic along the PE chain only."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import subprocess
import sys
import numpy as np
from ase.io import read, write

root = Path(__file__).resolve().parent
atoms = read(root / "model" / "pe_initial.extxyz")
out = root / "sampling" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
out.mkdir(parents=True)
env = {**os.environ, "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
results = []
for n in [1, 2, 4, 6]:
    case = out / f"repeat_c_{n}"
    case.mkdir()
    geometry = atoms.repeat((1, 1, n))
    write(case / "input.extxyz", geometry)
    with (case / "native_output.log").open("w") as log:
        subprocess.run([sys.executable, str(root / "run_gfn2.py"), "--worker", str(case / "input.extxyz"), str(case / "result.json")],
                       stdout=log, stderr=subprocess.STDOUT, check=True, env=env, timeout=600)
    data = json.loads((case / "result.json").read_text())
    force = np.asarray(data["forces_eV_per_A"]).reshape(n, len(atoms), 3)
    row = {"repeat_c": n, "atoms": len(geometry), "energy_eV": data["energy_eV"],
           "energy_per_primitive_cell_eV": data["energy_eV"]/n,
           "energy_per_CH2_eV": data["energy_eV"]/(4*n),
           "max_force_norm_eV_per_A": data["max_force_norm_eV_per_A"],
           "primitive_force_block_eV_per_A": force[0].tolist(),
           "max_force_replica_difference_eV_per_A": float(np.max(np.abs(force-force[0]))),
           "wall_seconds": data["wall_seconds"]}
    if results:
        row["energy_change_per_CH2_eV_from_previous"] = row["energy_per_CH2_eV"] - results[-1]["energy_per_CH2_eV"]
        row["max_force_change_from_previous_eV_per_A"] = float(np.max(np.abs(force[0]-np.asarray(results[-1]["primitive_force_block_eV_per_A"]))))
    results.append(row)
    summary = {"run_directory": str(out), "method": "GFN2-xTB via tblite 0.7 / ASE",
               "results": results, "scope": "Only c-direction supercell sampling at fixed repeated geometry. Transverse and production convergence untested.",
               "full_convergence_claim": False}
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    (root / "latest_sampling.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:v for k,v in row.items() if "block" not in k}), flush=True)
