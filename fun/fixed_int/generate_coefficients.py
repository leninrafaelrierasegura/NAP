"""Compute fixed-interval BRASIL coefficients; JSON request/result interface.

Setup: python3 -m venv --system-site-packages .venv-fixed-int
       .venv-fixed-int/bin/python -m pip install -r fun/fixed_int/requirements.txt
Probe: .venv-fixed-int/bin/python fun/fixed_int/generate_coefficients.py --check
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def toolchain():
    import numpy
    import scipy
    import baryrat
    import flamp
    import gmpy2
    from importlib.metadata import version
    return dict(python=sys.version.split()[0], executable=sys.executable,
                numpy=numpy.__version__, scipy=scipy.__version__,
                baryrat=baryrat.__version__, flamp=version("flamp"),
                gmpy2=gmpy2.version())


def fit_one(job, request):
    import numpy as np
    import baryrat
    import flamp
    from gmpy2 import mpfr
    digits = request["digits"]
    flamp.set_dps(digits)
    m, beta_text = job
    beta, lower = mpfr(beta_text), mpfr(request["x_min"])
    started = time.monotonic()
    target = lambda x: x**beta
    rational, info = baryrat.brasil(
        target, (lower, mpfr(1)), m, tol=request["tolerance"],
        maxiter=request["max_iterations"], npi=-40, info=True)
    if not info.converged:
        raise ValueError(f"BRASIL did not converge: m={m}, beta={beta_text}; "
                         "increase max_iterations or digits")
    # Check signed equioscillation as well as the algorithm's amplitude test.
    boundaries = np.concatenate(([lower], info.nodes, [mpfr(1)]))
    extrema, _ = baryrat.local_maxima_golden(
        lambda x: abs(target(x)-rational(x)), boundaries, num_iter=40)
    signed = target(extrema)-rational(extrema)
    normalized = signed / ((-1)**np.arange(len(signed))*np.sign(signed[0])*max(abs(signed)))
    equioscillation_gap = float(max(abs(normalized-1)))
    if equioscillation_gap > request["tolerance"]*2:
        raise ValueError(f"Signed equioscillation failed: m={m}, beta={beta_text}")
    # B(x) has poles t and residues v. S(y)=B(1/y) has poles 1/t,
    # residues -v/t^2 and constant B(0). Convert BEFORE rounding to double.
    tx, vx = rational.polres()
    imaginary = max(float(abs(x.imag)) for x in list(tx)+list(vx))
    if imaginary > 10**(-digits//2):
        raise ValueError(f"Nonreal source poles/residues: m={m}, beta={beta_text}")
    tx = np.array([mpfr(x.real) for x in tx], dtype=object)
    vx = np.array([mpfr(x.real) for x in vx], dtype=object)
    poles, residues = 1/tx, -vx/tx**2
    order = np.argsort(poles)
    poles, residues = poles[order], residues[order]
    constant = mpfr(rational(mpfr(0)))
    if (len(poles) != m or not all(poles < 0) or not all(residues > 0)
            or constant < 0 or len(set(poles)) != m):
        raise ValueError(f"Positive partial-fraction structure failed: m={m}, beta={beta_text}")
    xx = np.array([mpfr(0)]+[mpfr(str(x)) for x in
                  np.geomspace(float(lower), 1, 257)], dtype=object)
    bary_values = rational(xx)
    partial_mp = constant + np.array([
        sum(residues*x/(1-poles*x)) for x in xx], dtype=object)
    mp_gap = float(max(abs(partial_mp-bary_values)))
    pd, rd, kd, xd = np.array(poles, float), np.array(residues, float), float(constant), np.array(xx, float)
    double_values = kd+(xd[:, None]/(1-xd[:, None]*pd)).dot(rd)
    double_gap = float(max(abs(double_values-np.array(bary_values, float))))
    if mp_gap > 10**(-digits/2) or double_gap > 1e-12:
        raise ValueError(f"Partial-fraction conversion failed: m={m}, beta={beta_text}")
    return dict(m=m, beta=beta_text, x_min=request["x_min"],
        poles=pd.tolist(), residues=rd.tolist(), constant=kd,
        diagnostics=dict(converged=True, iterations=int(info.iterations),
            fit_max_abs_error=float(info.error), deviation=float(info.deviation),
            equioscillation_gap=equioscillation_gap, max_imaginary=imaginary,
            mp_conversion_gap=mp_gap, double_conversion_gap=double_gap,
            seconds=time.monotonic()-started),
        high_precision=dict(poles=[str(x) for x in poles], residues=[str(x) for x in residues],
            constant=str(constant), nodes=[str(x) for x in rational.nodes],
            values=[str(x) for x in rational.values], weights=[str(x) for x in rational.weights]))


def generate(request, workers):
    if request.get("schema_version") != 1:
        raise ValueError("Unsupported request schema")
    if not 0 < float(request["x_min"]) < 1:
        raise ValueError("x_min must lie strictly between 0 and 1")
    if not isinstance(request["digits"], int) or request["digits"] < 30:
        raise ValueError("digits must be an integer of at least 30")
    if (not 0 < request["tolerance"] < 1 or request["max_iterations"] < 1
            or workers < 1):
        raise ValueError("Invalid fitting controls")
    orders, betas = request["orders"], request["betas"]
    if any(not isinstance(m, int) or m < 1 for m in orders):
        raise ValueError("Orders must be positive integers")
    if any(not math.isfinite(float(b)) or not 0 < float(b) < 1 for b in betas):
        raise ValueError("Betas must lie strictly between 0 and 1")
    if len(set(orders)) != len(orders) or len(set(betas)) != len(betas):
        raise ValueError("Duplicate orders or betas")
    jobs = [(m, beta) for m in orders for beta in betas]
    fits = []
    if workers == 1:
        for job in jobs:
            fits.append(fit_one(job, request))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(fit_one, job, request): job for job in jobs}
            for count, future in enumerate(as_completed(futures), 1):
                try:
                    fits.append(future.result())
                except Exception:
                    for pending in futures:
                        pending.cancel()
                    raise
                if count % 10 == 0 or count == len(jobs):
                    print(f"BRASIL: {count}/{len(jobs)} fits completed", flush=True)
    fits.sort(key=lambda row: (float(row["beta"]), row["m"]))
    return dict(schema_version=1, request=request, toolchain=toolchain(), fits=fits)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--request")
    parser.add_argument("--output")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    try:
        if args.check:
            print(json.dumps(toolchain()))
            return
        if not args.request or not args.output:
            parser.error("--request and --output are required")
        result = generate(json.loads(Path(args.request).read_text()), args.workers)
        Path(args.output).write_text(json.dumps(result, allow_nan=False))
    except Exception as exc:
        print(f"BRASIL generation failed: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
