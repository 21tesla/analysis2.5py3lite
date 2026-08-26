"""
Command-line access to the NEF (BMRB NMR Enhanced Format, v1.1) import and
export implemented by ``ccpnmr/v2io/NefIo.py`` (read side) and
``ccpnmr/nefExport.py`` (write side).

Console entry point ``ccpnmr-nef`` (pyproject ``[project.scripts]``):

  ccpnmr-nef import <file.nef> [--project-name NAME] [--pdb PDB ...] [--force]
      Create a new CCPN project from the NEF file.  The project directory
      is created in the current working directory (named after the NEF
      file unless --project-name is given) and can then be opened in the
      GUI via Project > Open Project.  Optional PDB file(s) supply
      coordinates: one file reads all its models, several files read the
      first model of each (ensemble style).

  ccpnmr-nef export <project-directory> <output.nef>
      Write the NEF v1.1 file (metadata + molecular system + chemical
      shifts + restraints + peak lists) for an existing CCPN project.

NEF files carry metadata, peaks, shifts and restraints - never raw
spectrum matrix data.
"""
import argparse
import sys

__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2026"
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, https://doi.org/10.1007/s10858-016-0060-y")


def importNefFile(nefFilePath, projectName=None, pdbFilePaths=None, removeExisting=False):
    """Create a new CCPN project from a NEF file (NefIo.loadProject path).

    The project is saved to disk before returning, so it can be opened in
    the GUI (Project > Open Project) or re-processed by the export
    subcommand.  Returns the path of the created project directory (the
    'userData' repository).  Raises OSError if the project directory
    already exists and removeExisting is False.
    """
    from ccpnmr.v2io import NefIo
    from memops.general import Io as memopsIo

    memopsRoot = NefIo.loadProject(
        nefFilePath,
        pdbFilePaths=pdbFilePaths,
        projectName=projectName,
        removeExisting=removeExisting,
    )
    memopsIo.saveProject(memopsRoot)
    repo = memopsRoot.findFirstRepository(name="userData")
    return repo.url.path if repo is not None else None


def exportNefProject(projectDir, outFilePath):
    """Write the NEF (contemporary, v1.1) file for a CCPN project directory.

    Returns outFilePath.
    """
    from ccpnmr import nefExport
    from memops.general import Io as memopsIo

    memopsRoot = memopsIo.loadProject(projectDir)
    return nefExport.exportProject(memopsRoot, outFilePath)


def main(argv=None):
    """Console entry point (ccpnmr-nef). Returns the process exit code."""
    parser = argparse.ArgumentParser(
        prog="ccpnmr-nef",
        description="Import/export NEF (BMRB NMR Enhanced Format v1.1) project files",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_import = subparsers.add_parser("import", help="create a new CCPN project from a NEF file")
    p_import.add_argument("nef_file", help="path to the NEF file")
    p_import.add_argument(
        "--project-name",
        dest="project_name",
        default=None,
        help="project directory name (default: the NEF file name)",
    )
    p_import.add_argument(
        "--pdb",
        dest="pdb_files",
        nargs="+",
        default=None,
        metavar="PDB",
        help="PDB file(s) with coordinates to read into the project",
    )
    p_import.add_argument(
        "--force",
        action="store_true",
        help="delete the project directory if it already exists",
    )

    p_export = subparsers.add_parser("export", help="write a NEF file for a CCPN project directory")
    p_export.add_argument("project_dir", help="path to the project directory")
    p_export.add_argument("out_file", help="path of the NEF file to write")

    args = parser.parse_args(argv)

    # one --pdb file: read all models of that file; several: first model each
    pdbFilePaths = None
    if args.command == "import":
        if args.pdb_files:
            pdbFilePaths = args.pdb_files[0] if len(args.pdb_files) == 1 else args.pdb_files

    try:
        if args.command == "import":
            projectDir = importNefFile(
                args.nef_file,
                projectName=args.project_name,
                pdbFilePaths=pdbFilePaths,
                removeExisting=args.force,
            )
            print(f"wrote project directory {projectDir}")
            print("open it in the GUI via Project > Open Project")
        else:
            outFilePath = exportNefProject(args.project_dir, args.out_file)
            print(f"wrote {outFilePath}")
    except Exception as es:
        print(f"error: {es}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
