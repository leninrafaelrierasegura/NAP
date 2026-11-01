"""High-precision BRASIL fits for audit_recomputed_brasil.R.

Dependencies: numpy, scipy, gmpy2, baryrat==2.1.2, flamp==1.0.1.
No package installation is performed by this script. Optional --deps-dir adds
a directory of already installed dependencies to sys.path.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def fit_job(job, deps_dir, digits):
    if deps_dir:
        sys.path.insert(0, deps_dir)
    import numpy as np
    import baryrat
    import flamp
    import gmpy2
    from gmpy2 import mpfr
    flamp.set_dps(digits)
    started = time.monotonic()
    m = int(job["m"])
    beta = mpfr(job["beta"])
    lower = mpfr(job["x_min"])
    rational, info = baryrat.brasil(
        lambda x: x**beta, (lower, mpfr(1)), m,
        tol=1e-8, maxiter=2000, npi=-40, info=True)
    # If B(x) has poles t_j and residues v_j, then S(y)=B(1/y) has
    # p_j=1/t_j, r_j=-v_j/t_j^2 and k=B(0). Compute this conversion in MP,
    # and only then round its coefficients to doubles for the existing solver.
    tx, vx = rational.polres()
    max_imaginary = max(float(abs(x.imag)) for x in list(tx)+list(vx))
    if max_imaginary > 10**(-digits//2):
        raise ValueError(f"Complex poles/residues: {max_imaginary}")
    tx = np.array([mpfr(x.real) for x in tx], dtype=object)
    vx = np.array([mpfr(x.real) for x in vx], dtype=object)
    p = 1/tx
    residues = -vx/tx**2
    constant = mpfr(rational(mpfr(0)))
    order = np.argsort(p)
    p, residues = p[order], residues[order]
    valid = (len(p) == m and all(p < 0) and all(residues > 0)
             and constant >= 0 and len(set(p)) == m)
    # Conversion audit includes the fitting interval and x=0. The MP gap
    # distinguishes conversion error from approximation error.
    xx = np.concatenate(([mpfr(0)], np.array(
        [mpfr(str(x)) for x in np.geomspace(float(lower), 1, 257)], dtype=object)))
    bary_values = rational(xx)
    partial_mp = constant + np.array([
        sum(residues*x/(1-p*x)) for x in xx], dtype=object)
    conversion_mp_gap = float(max(abs(partial_mp-bary_values)))
    pd = np.array(p, dtype=float)
    rd = np.array(residues, dtype=float)
    kd = float(constant)
    xd = np.array(xx, dtype=float)
    partial_double = kd + (xd[:, None]/(1-xd[:, None]*pd)).dot(rd)
    double_gap = float(max(abs(partial_double-np.array(bary_values, dtype=float))))
    # Retain the barycentric representation and MP coefficients for inspection.
    serial = dict(job, digits=digits, converged=bool(info.converged),
        deviation=float(info.deviation), fit_max_abs_error=float(info.error),
        iterations=int(info.iterations), valid_positive_structure=bool(valid),
        max_imaginary=max_imaginary, mp_conversion_gap=conversion_mp_gap,
        double_conversion_gap=double_gap, seconds=time.monotonic()-started,
        baryrat_version=baryrat.__version__,
        nodes=[str(x) for x in rational.nodes], values=[str(x) for x in rational.values],
        weights=[str(x) for x in rational.weights],
        poles=[str(x) for x in p], residues=[str(x) for x in residues], constant=str(constant))
    coefficients = [dict(family=job["family"], beta=job["beta"], m=m, j=i+1,
        pole=float(p[i]), residue=float(residues[i]), constant=kd) for i in range(m)]
    return serial, coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--deps-dir", default="")
    parser.add_argument("--digits", type=int, default=70)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    with open(args.jobs, newline="") as stream:
        jobs = list(csv.DictReader(stream))
    results, coefficients = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fit_job, job, args.deps_dir, args.digits): job for job in jobs}
        for number, future in enumerate(as_completed(futures), 1):
            job = futures[future]
            try:
                result, rows = future.result()
            except Exception as exc:
                for pending in futures:
                    pending.cancel()
                raise RuntimeError(f"Fit failed for {job}: {exc}") from exc
            results.append(result)
            coefficients.extend(rows)
            # Checkpoint exact MP output as each job finishes.
            with (output/"fits.jsonl").open("a") as stream:
                stream.write(json.dumps(result)+"\n")
            if number % 10 == 0 or number == len(jobs):
                print(f"Completed {number}/{len(jobs)} fits; "
                      f"last={job['family']} beta={job['beta']} m={job['m']}", flush=True)
    coefficients.sort(key=lambda row: (row["family"], float(row["beta"]), row["m"], row["j"]))
    with (output/"coefficients.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(coefficients[0]))
        writer.writeheader()
        writer.writerows(coefficients)
    fields = [key for key in results[0] if key not in
              ("nodes", "values", "weights", "poles", "residues", "constant")]
    with (output/"fit_diagnostics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in results)


if __name__ == "__main__":
    main()
