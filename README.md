# Installing and Running CCPNMR Analysis 2.5.2 (Portable macOS)

This guide explains how to install and run the portable, self-contained standalone distribution of CCPNMR Analysis 2.5.2 on macOS (modernized for Python 3.13).

The standalone archive is fully self-contained. It embeds its own private Python 3.13 runtime, compiled C/Cython extension libraries, and all required scientific and graphical dependencies. **No developer tools, Xcode, or system python configurations are needed.**

---

## Prerequisites

CCPNMR Analysis uses a Tkinter/OpenGL interface, which requires **XQuartz** to run on macOS:

1. **Install XQuartz** (if you don't already have it):
   * **Via Homebrew**: `brew install --cask xquartz`
   * **Via Web Installer**: Download and install the `.pkg` from [xquartz.org](https://www.xquartz.org).
2. **Restart your Mac or Log out and back in** after the installation to initialize X11 services.

---

## Installation

1. Download the portable tarball (e.g., `ccpnmr-macos-arm64-standalone.tar.gz`).
2. Extract the archive in your Terminal:
   ```bash
   tar -xzf ccpnmr-macos-arm64-standalone.tar.gz
   cd ccpnmr-macos
   ```

---

## Running the Application

Use the provided `./run-ccpnmr.sh` launcher script. It will verify XQuartz is active and launch the respective application.

### 1. Launch the main CCPNMR Analysis GUI (Default)
```bash
./run-ccpnmr.sh
```

### 2. Launch other CCPN Utilities
The single launcher wrapper also serves all other sub-applications by passing a command argument:

* **Project Data Shifter:**
  ```bash
  ./run-ccpnmr.sh data-shifter
  ```
* **Project Format Converter:**
  ```bash
  ./run-ccpnmr.sh format-converter
  ```
* **Project Updater (Non-GUI):**
  ```bash
  ./run-ccpnmr.sh update
  ```

To see the help menu, run:
```bash
./run-ccpnmr.sh help
```

---

# Installing and Running CCPNMR Analysis 2.5.2 (Portable Linux x86_64)

The Linux standalone distribution is self-contained in the same way: an
embedded private CPython 3.13 runtime, the compiled C extensions, and all
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
