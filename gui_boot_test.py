#!/usr/bin/env python
"""P4-5 GUI boot test: launch each console-script app under Xvfb and assert it
reaches its mainloop (boots without raising) — the headless GUI gate.

Method: for each GUI app, run its `main()` entry in a subprocess under
`xvfb-run -a`. A GUI app that boots successfully enters `mainloop()` and keeps
running, so the run is killed after --timeout seconds and counted as BOOT OK.
An app that raises during import/construct exits early with a non-zero code and
a traceback — counted as FAILED (those are the real py3 runtime bugs this gate
exists to catch).

Usage:
    .venv/bin/python gui_boot_test.py                 # all apps, project venv
    .venv/bin/python gui_boot_test.py --apps ccpnmr,data-shifter
    /tmp/ccp-dist-venv/bin/python gui_boot_test.py    # installed-state run
Requires: xvfb-run on PATH (apt: xvfb / conda-forge: xvfb).
"""
import argparse
import os
import shutil
import subprocess
import sys

APPS = [
    # (name, call expression executed inside the app's python process)
    ("ccpnmr",
     "import ccpnmr.analysis.AnalysisGui as M; M.main()"),
    ("data-shifter",
     "import ccpnmr.format.gui.DataShifter as M; M.main([])"),
    ("format-converter",
     "import ccpnmr.format.gui.FormatConverter as M; M.main([])"),
]
NON_GUI = []

# SOURCE-STATE support: when this script sits next to the source tree, put that
# tree on the app subprocesses' PYTHONPATH (installed-state runs find the
# dist in site-packages and this is a harmless no-op).
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
APP_SUBPROCESS_ENV = dict(os.environ, MPLBACKEND="Agg")
_SRC_ROOT = os.path.join(_SCRIPT_DIR, "ccpnmr2.5", "python")
if os.path.isdir(_SRC_ROOT):
    _pp = APP_SUBPROCESS_ENV.get("PYTHONPATH", "")
    APP_SUBPROCESS_ENV["PYTHONPATH"] = _SRC_ROOT + (os.pathsep + _pp if _pp else "")


def run_app(name, code, python, timeout):
    if not code.startswith("import inspect"):
        cmd = ["xvfb-run", "-a", python, "-c", code]
    else:
        cmd = [python, "-c", code]
    try:
        proc = subprocess.run(
            cmd, env=APP_SUBPROCESS_ENV, timeout=timeout,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
    except subprocess.TimeoutExpired as e:
        # Still running at deadline => it reached mainloop / a dialog: boot OK.
        # (xvfb-run spawns Xvfb + child; kill the whole process group for hygiene.)
        killed_tail = (e.output or b"").decode("utf-8", "replace")[-400:]
        if "Traceback" in killed_tail:
            return False, "traceback while waiting (late startup crash)", killed_tail
        return True, "booted (killed at %.0fs in mainloop/dialog)" % timeout, killed_tail
    out = proc.stdout.decode("utf-8", "replace")
    if proc.returncode == 0 and "Traceback" not in out:
        return True, "exited cleanly (rc=0)", out[-400:]
    return False, "early exit rc=%s" % proc.returncode, out[-800:]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apps", default=None,
                    help="comma-separated subset of: %s (default: all)"
                         % ",".join(n for n, _ in APPS))
    ap.add_argument("--python", default=sys.executable,
                    help="python to run the apps with (default: this interpreter)")
    ap.add_argument("--timeout", type=float, default=15.0,
                    help="seconds to let a booted app run before killing (default 15)")
    args = ap.parse_args(argv)

    if shutil.which("xvfb-run") is None:
        print("ERROR: xvfb-run not found on PATH — install it (apt: xvfb / conda-forge: xvfb)")
        return 2

    selected = APPS if not args.apps else [a for a in APPS if a[0] in args.apps.split(",")]
    if not selected:
        print("ERROR: no apps matched --apps=%r" % args.apps)
        return 2

    print("GUI boot test (xvfb-run, timeout %.0fs, python=%s)" % (args.timeout, args.python))
    failures = []
    for name, code in selected:
        ok, why, tail = run_app(name, code, args.python, args.timeout)
        mark = "PASS" if ok else "FAIL"
        print("  [%s] %-16s — %s" % (mark, name, why))
        if not ok and "Traceback" in tail:
            tb = tail.split("Traceback", 1)[1]
            print("           " + " ".join(tb.split())[:400])
        if not ok:
            failures.append(name)
    # non-GUI entries: plain import/signature check (no xvfb needed)
    for name, code in NON_GUI:
        try:
            proc = subprocess.run([args.python, "-c", code], env=APP_SUBPROCESS_ENV,
                                  timeout=60, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT)
            ok = proc.returncode == 0 and b"Traceback" not in proc.stdout
            tail = proc.stdout.decode("utf-8", "replace")
        except subprocess.TimeoutExpired:
            ok, tail = False, "timeout"
        mark = "PASS" if ok else "FAIL"
        print("  [%s] %-40s — %s" % (mark, name, (tail.strip().splitlines() or ["?"])[-1][:120]))
        if not ok:
            failures.append(name)

    total = len(selected) + len(NON_GUI)
    print("%d/%d apps booted OK" % (total - len(failures), total))
    if failures:
        print("FAILED:", ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
