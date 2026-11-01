"""Run the supplied generator unchanged and audit its missing normalization.

Run from NAP with a Python containing numpy, scipy and baryrat:
  python test/brasil/audit_original_generators.py --deps-dir /tmp/nap-brasil-deps
The three supplied originals are preserved. Their CSVs are written into a
separate case directory, preventing collisions between beta values.
"""
import argparse
import contextlib
import csv
import hashlib
import io
import os
from pathlib import Path
import runpy
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deps-dir", default="")
    args = parser.parse_args()
    if args.deps_dir:
        sys.path.insert(0, args.deps_dir)
    import numpy as np
    import scipy.io
    import baryrat
    root = Path(__file__).resolve().parent
    output = root/"generator_audit"
    output.mkdir(exist_ok=True)
    package_file = output/"installed_coefficients.csv"
    package = list(csv.DictReader(package_file.open())) if package_file.exists() else []
    summaries = []
    for m, beta in [(8,.25), (8,.65), (8,.95), (14,.95), (15,.95), (20,.85)]:
        case = output/f"m{m}_beta{beta:.2f}"
        case.mkdir(exist_ok=True)
        old_cwd, old_args = Path.cwd(), sys.argv[:]
        log = io.StringIO()
        try:
            os.chdir(case)
            sys.argv = ["ratapprox.py", str(beta), str(m)]
            with contextlib.redirect_stdout(log):
                namespace = runpy.run_path(str(root/"ratapprox.py"), run_name="__main__")
        finally:
            os.chdir(old_cwd)
            sys.argv = old_args
        (case/"original_stdout.txt").write_text(log.getvalue())
        rational = namespace["r"]
        poles = np.loadtxt(case/f"poles{m}.csv", dtype=complex)
        zeros = np.loadtxt(case/f"zeros{m}.csv", dtype=complex)
        # fliplr(poly(zeros))/fliplr(poly(poles)) in the MATLAB file
        # evaluates at y=1 to the following ratio; its multiplicative gain is 1.
        legacy_at_one = np.prod(1-zeros)/np.prod(1-poles)
        gain = rational.gain()
        true_at_one = float(rational(1.0))
        y = np.geomspace(1, 10**((5+m)/2), 2001)
        actual = rational(1/y)
        legacy = np.prod(1-y[:,None]*zeros, axis=1)/np.prod(1-y[:,None]*poles, axis=1)
        subset = [row for row in package if int(row["m"]) == m and abs(float(row["beta"])-beta)<1e-12]
        installed_at_one = float("nan")
        if subset:
            installed_at_one = float(subset[0]["constant"])+sum(
                float(row["residue"])/(1-float(row["pole"])) for row in subset)
        row = dict(m=m, beta=beta, gain=float(np.real(gain)),
            original_barycentric_at_one=true_at_one,
            root_only_MATLAB_at_one=float(legacy_at_one.real),
            installed_table_at_one=installed_at_one,
            missing_gain_identity_gap=float(abs(legacy_at_one*gain-true_at_one)),
            gain_corrected_root_vs_barycentric_max_gap=float(max(abs(gain*legacy-actual))),
            maximum_root_imaginary=float(max(np.max(abs(poles.imag)), np.max(abs(zeros.imag)))),
            optimizer_printed_warning="warning:" in log.getvalue().lower())
        summaries.append(row)
        # MAT file permits checking MATLAB's exact poly()/residue() path without
        # starting the obsolete Mac-specific Python wrapper or the 999-row loop.
        scipy.io.savemat(case/"conversion_inputs.mat", dict(
            poles=poles, zeros=zeros, gain=gain, y=y, actual=actual, m=m, beta=beta))
        print(row, flush=True)
    with (output/"normalization.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    hashes = [f"{hashlib.sha256((root/name).read_bytes()).hexdigest()}  {name}"
              for name in ["ratapprox.py", "brasil.m", "rational_coeffs.m"]]
    (output/"original_sha256.txt").write_text("\n".join(hashes)+"\n")
    (output/"versions.txt").write_text(f"Python {sys.version}\nbaryrat {baryrat.__version__}\n")


if __name__ == "__main__":
    main()
