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

## NEF Project Files

CCPNMR Analysis reads and writes **NEF v1.1** (BMRB *Nmr_Exchange_Format*)
project files — metadata, molecules, chemical shifts, restraints and peak lists
(never raw spectrum matrix data):

* **GUI:** *Project → Load NEF…* creates a new CCPN project from a `.nef` file;
  *Project → Export NEF…* writes the current project to a `.nef` file.
* **Command line** (source / virtualenv installs, where the `ccpnmr-nef`
  console command is available):
  ```bash
  ccpnmr-nef import file.nef [--project-name NAME] [--pdb PDB ...] [--force]
  ccpnmr-nef export <project-directory> <output.nef>
  ```
