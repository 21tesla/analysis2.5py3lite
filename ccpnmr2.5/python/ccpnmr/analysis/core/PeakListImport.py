"""
======================COPYRIGHT/LICENSE START==========================

NmrdrawImport.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

Import of NMRdraw .tab peak list files into the NMR data model
(the "Add Peaks" feature).  The tabular format is the classic
VARS/FORMAT column format, e.g.::

    VARS   INDEX X_AXIS Y_AXIS ... X_PPM Y_PPM ... TYPE ASS CLUSTID MEMCNT
    FORMAT %5d %9.3f ...
    NULLVALUE -666
    NULLSTRING *
        1    27.858    17.204 ... 10.213 130.039 ... 1 W81-HE1    1    1

Each row is a peak (per-dimension PPM values in the ``<AXIS>_PPM``
columns).  The ASS column carries the atom assignment, e.g. ``W81-HE1``
(a 2 letter chain prefix, e.g. ``A-W81-HE1``, is also accepted).  For
every row a peak is added to a new peak list of the given data source,
and for every chemical shift a resonance is created and assigned to the
atom - the proton named in ASS for its own dimension, and the bonded
non-proton neighbour matching the dimension isotope for the other
dimensions (e.g. the 15N of a ``Q15-HN`` peak goes to atom N).  Rows that
cannot be mapped to an atom are ignored with a notification on stdout.
No second resonance is ever created for an (atom, isotope) that already
carries one - when that is the case (an "assignment already exists") the
new peak dimension is linked to the existing resonance only when
``overwrite`` is set, otherwise the existing assignment is kept and the
new dimension is left unassigned (both cases are notified on stdout).
"""

import os
import re

from ccpnmr.analysis.core import AssignmentBasic
from ccpnmr.analysis.core import PeakBasic
from ccpnmr.v2io import Constants

ASS_PATTERN = re.compile(
    r"^(?:(?P<chain>[A-Za-z0-9]{1,2})-)?(?P<res>[A-Za-z]+)(?P<num>\d+)-(?P<atom>.+)$"
)

_ISOTOPE_ELEMENT_PATTERN = re.compile(r"\d+([A-Z][a-z]?)")


def parseTabFile(filePath):
    """
    Parse an NMRdraw .tab peak list file (VARS/FORMAT column format).

    .. describe:: Input

    Word (file path)

    .. describe:: Output

    (List of Words, List of Dicts) - the VARS header names and one dict
    per data row, mapping header name to value.  NULLVALUE/NULLSTRING
    sentinels become None.

    Raises ValueError if the header is not recognisable.
    """
    with open(filePath, "r") as fp:
        lines = [line.rstrip("\n") for line in fp]

    varNames = None
    nullValue = None
    nullString = None
    rows = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if varNames is None:
            if tokens[0].upper() != "VARS":
                raise ValueError("Not a NMRdraw .tab file (missing VARS line): %s" % filePath)
            varNames = tokens[1:]
            continue
        if stripped.upper().startswith("FORMAT"):
            continue
        if stripped.upper().startswith("NULLVALUE"):
            nullValue = tokens[1] if len(tokens) > 1 else "-666"
            continue
        if stripped.upper().startswith("NULLSTRING"):
            nullString = tokens[1] if len(tokens) > 1 else "*"
            continue
        if len(tokens) < len(varNames):
            raise ValueError("Malformed data row (expected %d columns, got %d): %s" % (len(varNames), len(tokens), stripped[:60]))
        row = {}
        for name, token in zip(varNames, tokens):
            if nullString is not None and token == nullString:
                row[name] = None
            elif nullValue is not None and token == nullValue:
                row[name] = None
            else:
                row[name] = token
        rows.append(row)

    if varNames is None or not rows:
        raise ValueError("Not a NMRdraw .tab file (no VARS line or no data rows): %s" % filePath)

    return varNames, rows


def parseAssString(ass):
    """
    Parse an NMRdraw ASS assignment string, e.g. ``W81-HE1`` or
    ``A-W81-HE1``.

    .. describe:: Input

    Word (ASS value)

    .. describe:: Output

    Dict with 'chain' (or None), 'res', 'num' (int) and 'atom', or None
    if the string cannot be parsed.
    """
    if ass is None:
        return None
    ass = ass.strip()
    if not ass or ass.lower() in ("none", "unassigned"):
        return None
    match = ASS_PATTERN.match(ass)
    if match is None:
        return None
    return {
        "chain": match.group("chain"),
        "res": match.group("res"),
        "num": int(match.group("num")),
        "atom": match.group("atom"),
    }


def _isotopeElement(isotopeCode):
    """The element symbol implied by an isotope code, e.g. '15N' -> 'N'."""
    if not isotopeCode or isotopeCode == "unknown":
        return None
    match = _ISOTOPE_ELEMENT_PATTERN.match(isotopeCode)
    return match.group(1) if match else None


def _chemAtomNonHNeighbours(chemAtom):
    """List of the non-proton chemAtoms bonded to a given chemAtom."""
    try:
        bonds = list(chemAtom.chemBonds)
    except Exception:
        return []
    neighbours = []
    for bond in bonds:
        for other in bond.chemAtoms:
            if other is not chemAtom and not other.isDeleted and other.elementSymbol != "H":
                if other not in neighbours:
                    neighbours.append(other)
    return neighbours


# one-letter -> CCPN three-letter residue code (for ASS strings)
_ASS_RESIDUE_CODES = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}


def _assResidueCode(res):
    """The three-letter residue code implied by an ASS residue token."""
    if len(res) == 1:
        return _ASS_RESIDUE_CODES.get(res.upper())
    entry = Constants.residueName2chemCompId.get(res.upper())
    return entry[1] if entry else None


def findTargetAtom(atom, isotopeCode):
    """
    Find the atom to assign a resonance of the given dimension isotope to,
    given the atom named in the ASS string of the peak.  The named atom is
    used when its element matches the isotope (the proton dimension);
    otherwise the bonded non-proton neighbour with a matching element is
    used (e.g. the 15N dimension of a HN peak goes to atom N).

    .. describe:: Input

    MolSystem.Atom, Word (isotope code, e.g. '1H')

    .. describe:: Output

    MolSystem.Atom or None
    """
    targetElement = _isotopeElement(isotopeCode)
    if targetElement is None:
        return None

    chemAtom = atom.chemAtom
    if chemAtom is not None and chemAtom.elementSymbol == targetElement:
        return atom

    for neighbour in _chemAtomNonHNeighbours(chemAtom) if chemAtom else []:
        if neighbour.elementSymbol == targetElement:
            return atom.residue.findFirstAtom(name=neighbour.name)

    return None


def _findResidue(molSystems, assInfo):
    """
    Find the residue matching a parsed ASS assignment.

    .. describe:: Output

    MolSystem.Residue or None (first match wins when several chains or
    mol systems carry the same seqCode)
    """
    expectedCode = _assResidueCode(assInfo["res"])
    for molSystem in molSystems:
        if molSystem.isDeleted:
            continue
        for chain in molSystem.sortedChains():
            if chain.isDeleted:
                continue
            if assInfo["chain"] is not None and chain.code != assInfo["chain"]:
                continue
            residue = chain.findFirstResidue(seqCode=assInfo["num"])
            if residue is None:
                residue = chain.findFirstResidue(seqCode=str(assInfo["num"]))
            if residue is None or residue.isDeleted:
                continue
            if expectedCode is not None and residue.ccpCode != expectedCode:
                continue
            return residue
    return None


def _resolveAtom(residue, atomName):
    """
    Resolve an ASS atom name to the residue's atom.  Besides the exact
    name match this understands the backbone amide proton: NMRdraw and
    the CCPN nomenclature call it HN, while some molecules (like the
    NMRdraw-derived ones) name the backbone H bonded to N simply "H".
    """
    atom = residue.findFirstAtom(name=atomName)
    if atom is not None:
        return atom
    if atomName == "HN":
        candidate = residue.findFirstAtom(name="H")
        if candidate is not None and candidate.chemAtom is not None:
            for neighbour in _chemAtomNonHNeighbours(candidate.chemAtom):
                if neighbour.elementSymbol == "N":
                    return candidate
    return None


def _resonanceGroupFor(residue, nmrProject, groups):
    """The per-residue ResonanceGroup, created on demand (NEF-importer style)."""
    group = groups.get(id(residue))
    if group is None:
        chain = residue.chain
        name = "%s.%s" % (chain.code, residue.seqCode)
        group = nmrProject.newResonanceGroup(name=name)
        group.molType = residue.molType
        group.ccpCode = residue.ccpCode
        group.residue = residue
        groups[id(residue)] = group
    return group


def _existingAtomResonance(atom, isotopeCode):
    """
    The resonance, if any, already assigned to this atom for this isotope,
    whether or not it is already linked to peak dimensions.  The import
    never creates a second resonance for the same (atom, isotope) - the
    Overwrite resonance option decides only whether the new peak dimension
    is additionally linked to it.

    .. describe:: Input

    MolSystem.Atom, Word (isotope code)

    .. describe:: Output

    Nmr.Resonance or None
    """
    atomSet = atom.atomSet
    if atomSet is None:
        return None
    for resonanceSet in atomSet.resonanceSets:
        if resonanceSet.isDeleted:
            continue
        for resonance in resonanceSet.resonances:
            if not resonance.isDeleted and resonance.isotopeCode == isotopeCode:
                return resonance
    return None


def _rowIntensity(row, column):
    """The row value of a numeric column as float, or None if absent/unusable."""
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError):
        return None
    if value != value or value in (float("inf"), float("-inf")):
        return None
    return value


def _makePeak(peakList, row, hasHeight, hasVolume, annotation=None):
    """Add a new peak to the list, with the row's intensity values if present."""
    parameters = {}
    if annotation:
        parameters["annotation"] = annotation
    height = _rowIntensity(row, "HEIGHT") if hasHeight else None
    volume = _rowIntensity(row, "VOL") if hasVolume else None
    if height is not None:
        parameters["height"] = height
    if volume is not None:
        parameters["volume"] = volume
    peak = peakList.newPeak(**parameters)
    # the helper functions must set the attributes - setting them directly
    # crashes in c-code (same requirement as the NEF importer)
    if height is not None:
        PeakBasic.setManualPeakIntensity(peak, height, intensityType="height")
    if volume is not None:
        PeakBasic.setManualPeakIntensity(peak, volume, intensityType="volume")
    return peak


def importTabPeaks(api, filePath, dataSource, listName=None, shiftList=None, overwrite=False):
    """
    Import the peaks of an NMRdraw .tab file or a .nef file into a new
    peak list of a data source, creating and assigning resonances.

    .. describe:: Input

    ApiBase (project loaded in memory), Word (.tab or .nef file path),
    Nmr.DataSource (spectrum to take the peaks + make the peak list for),
    Word (peak list name, default root of file name),
    Nmr.ShiftList (resonance list to populate, default the experiment's),
    Bool (overwrite existing resonance assignments)

    .. describe:: Output

    Dict report: file, peakList, rows, peaksAdded, peaksUnassigned,
    peaksSkipped (rows that could not be mapped), resonancesCreated,
    assignmentsApplied, assignmentsOverwritten, assignmentsKept,
    dimResonancesUnassigned (shifts with no mappable atom), error (None
    or Word).  Progress notifications are sent to stdout.
    """
    report = {
        "file": filePath,
        "peakList": None,
        "rows": 0,
        "peaksAdded": 0,
        "peaksUnassigned": 0,
        "peaksSkipped": 0,
        "resonancesCreated": 0,
        "assignmentsApplied": 0,
        "assignmentsOverwritten": 0,
        "assignmentsKept": 0,
        "dimResonancesUnassigned": 0,
        "error": None,
    }

    prefix = "Add Peaks: "
    isNef = filePath.lower().endswith(".nef")
    try:
        if isNef:
            from ccpnmr.nef import StarIo
            nefDataExtent = StarIo.parseNefFile(filePath)
            data_blocks = list(nefDataExtent.values())
            if not data_blocks:
                raise ValueError("No data blocks found in NEF file")
            data_block = data_blocks[0]
            
            spectrum_saveframes = [
                sf for sf in data_block.values()
                if isinstance(sf, StarIo.NmrSaveFrame) and sf.category == "nef_nmr_spectrum"
            ]
            if not spectrum_saveframes:
                raise ValueError("No spectrum saveframe found in NEF file")
            
            peak_saveframe = None
            ds_name = dataSource.name.lower() if dataSource.name else ""
            exp_name = dataSource.experiment.name.lower() if dataSource.experiment.name else ""
            
            for sf in spectrum_saveframes:
                sf_name = sf.name.lower()
                if (ds_name and ds_name in sf_name) or (exp_name and exp_name in sf_name):
                    peak_saveframe = sf
                    break
            if peak_saveframe is None:
                peak_saveframe = spectrum_saveframes[0]
                
            if "nef_peak" not in peak_saveframe:
                raise ValueError("No nef_peak loop found in spectrum saveframe")
                
            loop = peak_saveframe["nef_peak"]
            varNames = list(loop.columns)
            rows = loop.data
            
            # Map NEF dimension IDs to normalized isotope/axis codes
            nef_dims = {}
            if "nef_spectrum_dimension" in peak_saveframe:
                dim_loop = peak_saveframe["nef_spectrum_dimension"]
                for dim_row in dim_loop.data:
                    nef_dims[dim_row["axis_code"].strip().upper()] = int(dim_row["dimension_id"])
        else:
            varNames, rows = parseTabFile(filePath)
    except Exception as e:
        report["error"] = str(e)
        print(prefix + "ERROR %s" % e)
        return report

    if listName is None or not listName.strip():
        listName = os.path.splitext(os.path.basename(filePath))[0]
    if shiftList is None:
        shiftList = dataSource.experiment.shiftList
    if shiftList is None:
        shiftList = AssignmentBasic.getShiftLists(api.currentNmrProject)[0] if AssignmentBasic.getShiftLists(api.currentNmrProject) else None
    if shiftList is None:
        report["error"] = "No resonance (shift) list available - create one first"
        print(prefix + "ERROR no resonance (shift) list available")
        return report

    nmrProject = api.currentNmrProject
    if shiftList is not dataSource.experiment.shiftList:
        dataSource.experiment.shiftList = shiftList

    molSystems = [ms for ms in api.sortedMolSystems() if not ms.isDeleted]
    if not molSystems:
        report["error"] = "No molecule defined - use the Add Molecule button first"
        print(prefix + "ERROR no molecule defined")
        return report

    axisLetters = "XYZWUV"
    if isNef:
        missing = [dim for dim in range(1, dataSource.numDim + 1) if "position_%d" % dim not in varNames]
        if missing:
            report["error"] = "File is missing the position_%s column(s) for a %dD spectrum" % (", ".join(str(dim) for dim in missing), dataSource.numDim)
            print(prefix + "ERROR " + report["error"])
            return report
    else:
        missing = [letter for letter in axisLetters[: dataSource.numDim] if "%s_PPM" % letter not in varNames]
        if missing:
            report["error"] = "File is missing the %s_PPM column(s) for a %dD spectrum" % (", ".join("%s_PPM" % letter for letter in missing), dataSource.numDim)
            print(prefix + "ERROR " + report["error"])
            return report

    dataDims = dataSource.sortedDataDims()
    dataDimRefs = {dataDim.dim: dataDim.findFirstDataDimRef() for dataDim in dataDims}
    dimIsotopes = {}
    for dataDim in dataDims:
        ref = dataDimRefs[dataDim.dim]
        if ref is not None and ref.expDimRef is not None and ref.expDimRef.isotopeCodes:
            dimIsotopes[dataDim.dim] = ref.expDimRef.isotopeCodes[0]
        else:
            dimIsotopes[dataDim.dim] = "unknown"

    peakList = dataSource.newPeakList(name=listName, details=listName)
    report["peakList"] = peakList
    groups = {}
    seenPositions = {}
    hasHeight = "HEIGHT" in varNames or "height" in varNames
    hasVolume = "VOL" in varNames or "volume" in varNames

    print(prefix + "Importing %d line(s) of %s into peak list '%s'" % (len(rows), os.path.basename(filePath), listName))

    for row in rows:
        report["rows"] += 1
        rowNumber = row.get("INDEX") or row.get("index") or row.get("peak_id")

        if isNef:
            if "height" in row:
                row["HEIGHT"] = row["height"]
            if "volume" in row:
                row["VOL"] = row["volume"]
            try:
                positions = []
                for dataDim in dataDims:
                    iso = dimIsotopes[dataDim.dim].strip().upper()
                    nef_dim = nef_dims.get(iso, dataDim.dim)
                    positions.append(float(row["position_%d" % nef_dim]))
            except (KeyError, TypeError, ValueError):
                report["peaksSkipped"] += 1
                print(prefix + "row %s: unmappable (missing position values) - ignored" % rowNumber)
                continue
        else:
            try:
                positions = [float(row["%s_PPM" % axisLetters[dataDim.dim - 1]]) for dataDim in dataDims]
            except (KeyError, TypeError, ValueError):
                report["peaksSkipped"] += 1
                print(prefix + "row %s: unmappable (missing PPM values) - ignored" % rowNumber)
                continue

        key = tuple(round(value, 3) for value in positions)
        if key in seenPositions:
            report["peaksSkipped"] += 1
            print(prefix + "row %s: duplicate position %s - ignored" % (rowNumber, tuple(positions)))
            continue
        seenPositions[key] = True

        if isNef:
            has_any_assignment = False
            dim_assignments = {}
            dim_groups = {}
            annotation_parts = []
            
            for dataDim in dataDims:
                iso = dimIsotopes[dataDim.dim].strip().upper()
                nef_dim = nef_dims.get(iso, dataDim.dim)
                
                chain = row.get("chain_code_%d" % nef_dim)
                seq_code = row.get("sequence_code_%d" % nef_dim)
                res_name = row.get("residue_name_%d" % nef_dim)
                atom_name = row.get("atom_name_%d" % nef_dim)
                
                if seq_code in (None, ".", "?") or atom_name in (None, ".", "?"):
                    continue
                
                try:
                    seq_num = int(seq_code)
                except (TypeError, ValueError):
                    seq_num = seq_code
                
                dimAssInfo = {
                    "chain": chain if chain not in (".", "?") else None,
                    "num": seq_num,
                    "res": res_name if res_name not in (".", "?") else None,
                    "atom": atom_name,
                }
                
                residue = _findResidue(molSystems, dimAssInfo)
                atom = None
                if residue is not None:
                    atom = _resolveAtom(residue, atom_name)
                    
                if atom is not None:
                    has_any_assignment = True
                    dim_assignments[dataDim.dim] = atom
                    dim_groups[dataDim.dim] = _resonanceGroupFor(residue, nmrProject, groups)
                    annotation_parts.append("%s%s-%s" % (res_name, seq_code, atom_name))
            
            annotation = " ".join(annotation_parts) if annotation_parts else None
            
            if not has_any_assignment:
                peak = _makePeak(peakList, row, hasHeight, hasVolume)
                report["peaksUnassigned"] += 1
                print(prefix + "row %s: no assignment - peak at %s added unassigned" % (rowNumber, tuple(positions)))
                continue
                
            peak = _makePeak(peakList, row, hasHeight, hasVolume, annotation=annotation)
        else:
            ass = row.get("ASS")
            assInfo = parseAssString(ass)

            if assInfo is None:
                peak = _makePeak(peakList, row, hasHeight, hasVolume)
                report["peaksUnassigned"] += 1
                print(prefix + "row %s: no assignment (%s) - peak at %s added unassigned" % (rowNumber, ass, tuple(positions)))
                continue
            else:
                residue = _findResidue(molSystems, assInfo)
                atom = None
                if residue is not None:
                    atom = _resolveAtom(residue, assInfo["atom"])
                if atom is None:
                    report["peaksSkipped"] += 1
                    print(prefix + "row %s: unmappable (no atom %s on residue %s%d in the molecule) - ignored"
                          % (rowNumber, assInfo["atom"], assInfo["res"], assInfo["num"]))
                    continue

                peak = _makePeak(peakList, row, hasHeight, hasVolume, annotation=ass.strip())
                group = _resonanceGroupFor(residue, nmrProject, groups)

        for dataDim, position in zip(dataDims, positions):
            peakDim = None
            for candidate in peak.sortedPeakDims():
                if candidate.dim == dataDim.dim:
                    peakDim = candidate
                    break
            if peakDim is None:
                continue
            # the value setter needs an unambiguous dataDimRef first
            peakDim.dataDimRef = dataDimRefs[dataDim.dim]
            peakDim.value = position
            isotopeCode = dimIsotopes[dataDim.dim]
            
            if isNef:
                targetAtom = dim_assignments.get(dataDim.dim)
                group = dim_groups.get(dataDim.dim)
            else:
                targetAtom = findTargetAtom(atom, isotopeCode) if atom is not None else None

            if targetAtom is None:
                # A resonance is still created for the chemical shift,
                # but no atom can stand behind it for this dimension
                resonance = nmrProject.newResonance(
                    name=_isotopeElement(isotopeCode) or "unknown",
                    isotopeCode=isotopeCode,
                    resonanceGroup=group,
                )
                report["resonancesCreated"] += 1
                report["dimResonancesUnassigned"] += 1
                if isNef:
                    print(prefix + "row %s: dim %d has no atom assigned (position %s) - resonance left unassigned to an atom"
                          % (rowNumber, dataDim.dim, position))
                else:
                    print(prefix + "row %s: no %s atom behind ASS %s (position %s) - resonance left unassigned to an atom"
                          % (rowNumber, isotopeCode, ass, position))
                AssignmentBasic.newPeakDimContrib(peakDim, resonance)
                AssignmentBasic.updateResonShift(resonance, peakDim)
                continue

            existing = _existingAtomResonance(targetAtom, isotopeCode)
            if existing is not None:
                if not overwrite:
                    report["assignmentsKept"] += 1
                    if isNef:
                        print(prefix + "row %s: %s of %s already has a resonance (name %s, serial %s) - kept, new dim left unassigned (Overwrite resonance off)"
                              % (rowNumber, isotopeCode, targetAtom.name, existing.name, existing.serial))
                    else:
                        print(prefix + "row %s: %s of %s already has a resonance (name %s, serial %s) - kept, new dim left unassigned (Overwrite resonance off)"
                              % (rowNumber, isotopeCode, ass, existing.name, existing.serial))
                    continue
                if isNef:
                    print(prefix + "row %s: overwriting - reusing existing %s resonance (name %s, serial %s) for %s"
                          % (rowNumber, isotopeCode, existing.name, existing.serial, targetAtom.name))
                else:
                    print(prefix + "row %s: overwriting - reusing existing %s resonance (name %s, serial %s) for %s"
                          % (rowNumber, isotopeCode, existing.name, existing.serial, ass))
                resonance = existing
                report["assignmentsOverwritten"] += 1
            else:
                resonance = nmrProject.newResonance(
                    name=targetAtom.name,
                    isotopeCode=isotopeCode,
                    resonanceGroup=group,
                )
                report["resonancesCreated"] += 1
                atomSet = targetAtom.atomSet
                if atomSet is None:
                    atomSet = nmrProject.newAtomSet(atoms=[targetAtom])
                AssignmentBasic.assignAtomsToRes([atomSet], resonance)

            AssignmentBasic.newPeakDimContrib(peakDim, resonance)
            AssignmentBasic.updateResonShift(resonance, peakDim)
            report["assignmentsApplied"] += 1
        report["peaksAdded"] += 1

    print(prefix + "Done: %d peak(s) added (%d unassigned, %d row(s) skipped), %d resonance(s) created, "
          "%d assignment(s) applied (%d overwritten, %d kept)"
          % (report["peaksAdded"], report["peaksUnassigned"], report["peaksSkipped"],
             report["resonancesCreated"], report["assignmentsApplied"],
             report["assignmentsOverwritten"], report["assignmentsKept"]))

    return report
