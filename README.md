# Installing and Running CCPNMR Analysis 2.5.2 (Portable macOS)

This guide explains how to install and run the portable, self-contained standalone distribution of CCPNMR Analysis 2.5.2 on macOS (modernized for Python 3.14).

The standalone archive is fully self-contained. It embeds its own private Python 3.14 runtime, compiled C/Cython extension libraries, and all required scientific and graphical dependencies. **No developer tools, Xcode, or system python configurations are needed.**

---

## Prerequisites

**None.** The macOS standalone is fully self-contained: it embeds its own
private CPython 3.14 runtime together with Tcl/Tk 9.0.4 (native Aqua canvas
port — no X11 / XQuartz required). Just download, extract, and run.

---

## Installation

1. Download the portable tarball (e.g., `ccpnmr-2.5.2-macos-arm64-standalone.tar.gz`).
2. Extract the archive in your Terminal:
   ```bash
   tar -xzf ccpnmr-2.5.2-macos-arm64-standalone.tar.gz
   cd ccpnmr-2.5.2-macos-arm64-standalone
   ```

> The archive is built for the architecture of the Mac that built it
> (Apple Silicon = `arm64`, Intel = `x86_64`) — run it on a matching Mac.

---

## Running the Application

There is no launcher script — the tree is self-contained. Run:

### 1. Launch the main CCPNMR Analysis GUI (Default)
```bash
./bin/analysis                # or: ./bin/analysis /path/to/project
```

### 2. NEF import / export (non-GUI)
```bash
./runtime/bin/ccpnmr-nef import file.nef [--project-name NAME] [--force]
./runtime/bin/ccpnmr-nef export <project-directory> <output.nef>
```

The tree reads nothing from and writes nothing to the system except the project
you open/save (and your default browser for Project > Summary).

---

# Installing and Running CCPNMR Analysis 2.5.2 (Portable Linux x86_64)

The Linux standalone distribution is self-contained in the same way: an
embedded private CPython 3.14 runtime, the compiled C extensions, and all
dependencies. No system Python and no `pip` are involved at run time.

1. Unpack the standalone archive (e.g. `ccpnmr-2.5.2-linux-x86_64-standalone.tar.gz`):
   ```bash
   tar -xzf ccpnmr-2.5.2-linux-x86_64-standalone.tar.gz
   cd ccpnmr-2.5.2-linux-x86_64-standalone
   ```
2. Run it (optionally pointing at an existing project directory):
   ```bash
   ./bin/analysis                 # or: ./bin/analysis /path/to/project
   ```
3. Non-GUI utilities (NEF import/export):
   ```bash
   ./runtime/bin/ccpnmr-nef import file.nef [--project-name NAME] [--force]
   ./runtime/bin/ccpnmr-nef export <project-directory> <output.nef>
   ```

Host requirements: a Linux x86_64 desktop with the usual X11 graphics
libraries (`libX11`, `libGL`) — nothing else is read from or written to
the system. The tree is relocatable: move it anywhere after unpacking.

> The standalone tree is produced by `./make-standalone-linux.sh` in the
> source repository (it rebuilds the wheel and packs it with the private
> runtime into `dist/`).

---

## NEF Project Files

CCPNMR Analysis reads and writes **NEF v1.1** (BMRB *Nmr_Exchange_Format*)
project files — metadata, molecules, chemical shifts, restraints and peak lists
(never raw spectrum matrix data):

* **GUI:** *Project → Export NEF…* writes the current project to a `.nef`
  file (NEF is metadata + model, never raw spectrum matrix data — import is
  command-line only).
* **Command line** (standalone: `./runtime/bin/ccpnmr-nef`; source /
  virtualenv installs: the `ccpnmr-nef` console command):
  ```bash
  ccpnmr-nef import file.nef [--project-name NAME] [--pdb PDB ...] [--force] [--relink [DIR]]
  ccpnmr-nef export <project-directory> <output.nef>
  ```
