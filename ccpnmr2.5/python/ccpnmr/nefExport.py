"""Export a legacy (data model 2.1.2) NMR project to a contemporary NEF file.

Stage 37 write-side companion to ``ccpnmr/v2io/NefIo.py`` (the read-side
importer).  The pair covers reading AND writing contemporary NEF files
(BMRB NEF v1.1, ``Nmr_Exchange_Format``) against the legacy model; legacy
file handling is out of scope.

Contents exported
-----------------
saveframes are written in the reading order used by
``v2io.NefIo.saveFrameReadingOrder`` so that a file exported here imports
cleanly with ``v2io.NefIo.loadNefFile`` / ``loadProject``:

  nef_nmr_meta_data                  (mandatory)
  nef_molecular_system               (mandatory; ``nef_sequence`` loop)
  nef_chemical_shift_list_<name>     (one per ShiftList)
  nef_distance_restraint_list_<name> (HBondConstraintList → ``restraint_origin=hbond``)
  nef_dihedral_restraint_list_<name>
  nef_rdc_restraint_list_<name>
  ccpn_restraint_list_<name>         (JCoupling / ChemicalShift / Csa)
  nef_nmr_spectrum_<name>            (one per DataSource: dimensions,
                                      transfers, peak list + peaks,
                                      incl. per-dimension assignments)
  nef_peak_restraint_links           (only if constraints carry peaks)

Conventions
-----------
  * Column sets mirror a real CCPN-exported NEF v1.1 file (see the bundled
    ``ccpnmr/nef/testdata/CCPN_Commented_Example.nef``); ``.`` is written
    for undefined values (the parser reads this back as ``None``).
  * Resonance identity (chain_code / sequence_code / residue_name /
    atom_name) is derived from the Resonance + its ResonanceGroup, the
    inverse of ``NefIo.CcpnNefReader.fetchResidueMap``: an RG name
    ``'A.63'`` → chain ``A`` seq ``63``; an RG without a chain prefix is
    written under chain ``'@'`` (the NEF convention for unassigned
    chains).
  * A peak with N alternative assignments (N PeakContribs) is written as
    N rows sharing the same ``peak_id``, exactly the multi-row form the
    importer understands.
  * Raw spectrum matrix data is not part of NEF.  A DataSource linked
    to a data file additionally carries the reference (file path +
    file-format items + the data-dimension point counts) so a plain
    same-machine reimport auto-links it; unlinked DataSources export
    with no file fields (the NEF carries peaks + shifts + restraints
    + metadata either way).

Public API
----------
  makeNefDataBlock(memopsRoot) -> ccpnmr.nef.StarIo.NmrDataBlock
  exportProject(memopsRoot, fileName) -> fileName
"""

import os
import random
import sys

from ccpnmr import Common as commonUtil
from ccpnmr.nef import StarIo

__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2026"
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, https://doi.org/10.1007/s10858-016-0060-y")

programName = 'CcpNmr'
programVersion = '2.5.2'


# ---------------------------------------------------------------------------
# identity extraction (resonance, isotope, sequence)
# ---------------------------------------------------------------------------

def _residueName(residue):
    """Standard 3-letter residue_name for a MolResidue.

    The importer resolves residue_name through ``v2io.Constants.
    residueName2chemCompId``, whose keys are the standard codes (ALA, GLY,
    ...; DA/DC/DG/DT and A/G/C/U for nucleic acids); the legacy title-case
    ``ccpCode`` (Ala, Gly, ...; T) is NOT a valid key.  The standard code is
    therefore taken from the residue's ChemComp (``code3Letter``), the same
    form a real CCPN-exported NEF file carries."""
    chemCompVar = getattr(residue, 'chemCompVar', None)
    chemComp = getattr(chemCompVar, 'chemComp', None) if chemCompVar is not None else None
    if chemComp is not None:
        code = getattr(chemComp, 'code3Letter', None) or getattr(chemComp, 'ccpCode', None)
        if code:
            return code
    return 'UNK'


def _canonicalAtomName(name, mappingType):
    """NEF-canonical atom_name for a model AtomSetMapping name.

    Byte-for-byte the transform the importer applies when it builds its
    per-residue atomMappings keys (``NefIo.load_nef_sequence``): pseudoatom
    ``*`` -> ``%``, upper cased, and nonstereo prochiral endings A/B ->
    x/y (%%-suffixed first, since they never end in a bare letter).  A
    row carrying this name resolves to the mapped AtomSetMapping on
    reimport; the raw model name would not."""
    atName = name.replace("*", "%").upper()
    if mappingType == "nonstereo":
        for tag, val in (("A", "x"), ("B", "y"), ("A%", "x%"), ("B%", "y%")):
            if atName.endswith(tag):
                return atName[: -len(tag)] + val
    return atName


def _resonanceAtomNames(resonance):
    """The set of MolSystem atom names covered by a resonance's Set.

    None when the resonance carries no atom identity (no resonanceSet or
    no atom names on its sets)."""
    resonanceSet = getattr(resonance, 'resonanceSet', None)
    if resonanceSet is None:
        return None
    names = set()
    for atomSet in resonanceSet.atomSets:
        for atom in atomSet.atoms:
            atomName = getattr(atom, 'name', None)
            if atomName:
                names.add(atomName)
    return names or None


def _resonanceAtomName(resonance, residue):
    """Recover the NEF-canonical atom_name from a nameless resonance's atoms.

    In a natively created legacy project ``resonance.name`` is None - the
    atom identity lives only in the resonance's resonanceSet (atomSets of
    MolSystem atoms).  Matching those atoms against the residue's
    ResidueMapping atomSetMappings (the same mapping the importer builds
    and its ``fetchAtomMap`` keys on) yields the model name, which
    ``_canonicalAtomName`` then converts to the form the importer
    resolves.  Returns None when the resonance carries no usable atom
    identity (no resonanceSet, no matching mapping) so the caller can fall
    back to the element@serial pin."""
    if residue is None:
        return None
    resoAtoms = _resonanceAtomNames(resonance)
    if not resoAtoms:
        return None
    try:
        residueMapping = getattr(residue, 'residueMapping', None)
        if residueMapping is None:
            from ccpnmr.analysis.core import MoleculeBasic

            residueMapping = MoleculeBasic.getResidueMapping(residue, aromaticsEquivalent=True)
    except Exception:
        return None
    candidates = []
    for atomSetMapping in residueMapping.atomSetMappings:
        mappedAtoms = set()
        for atomSet in atomSetMapping.atomSets:
            for atom in atomSet.atoms:
                atomName = getattr(atom, 'name', None)
                if atomName:
                    mappedAtoms.add(atomName)
        if mappedAtoms and mappedAtoms == resoAtoms:
            name = getattr(atomSetMapping, 'name', None)
            if name:
                candidates.append(atomSetMapping)
    if not candidates:
        return None
    # one atom group can carry several mappings - e.g. the HE2/HE3 pair
    # has stereo (He2/He3), nonstereo (Hda/Hdb covering the pair) AND
    # ambiguous (He*) spellings, all resolving from the same atom names
    # (the ResidueMapping reuses ONE AtomSet across differently-named
    # asms, so neither object nor name-set identity can pick one); the
    # spelling choice must be deterministic AND import-safe, because one
    # NEF row per DISTINCT atom_name resolves one-for-one, a '%' row over
    # an ambiguous mapping EXPANDS on reimport (one row -> two resonances,
    # duplicating the row's shift), and a REPEATED row identity collides
    # (the second row resolves to the first's resonance, which already
    # holds a Shift in the list - the model allows one Shift per
    # (list, resonance)).  Hence:
    # - a LONE resonance over its group takes the nonstereo spelling
    #   (reimports one-for-one; writing it '%' would expand into a
    #   phantom second resonance), else the ambiguous '%' (which recreates
    #   exactly the pair a single real NEF '%' row carries);
    # - a SIBLING class (several nameless resonances over the same atom
    #   names, e.g. the two protons of a GLY Halpha CH2) gets DISTINCT
    #   import-resolvable spellings, one per sibling - nonstereo first
    #   (its combined atom sets recreate each sibling's own resonanceSet
    #   shape), then per-atom stereo; deterministic by resonance serial.
    #   A group offering no distinct spelling past the shared '%' form is
    #   left to the element@serial pin (import-safe).
    def _byName(asm):
        return (getattr(asm, 'name', '') or '')

    nonstereo = sorted(
        (asm for asm in candidates if getattr(asm, 'mappingType', None) == 'nonstereo'),
        key=_byName)
    stereo = sorted(
        (asm for asm in candidates if getattr(asm, 'mappingType', None) == 'stereo'),
        key=_byName)
    ambiguous = [asm for asm in candidates if getattr(asm, 'mappingType', None) == 'ambiguous']

    siblings = []
    group = getattr(resonance, 'resonanceGroup', None)
    for other in (list(getattr(group, 'resonances', None) or [])):
        if other is resonance:
            continue
        if getattr(other, 'isotopeCode', None) != getattr(resonance, 'isotopeCode', None):
            continue
        if _resonanceAtomNames(other) == resoAtoms:
            siblings.append(other)

    if siblings:
        ordered = sorted([resonance] + siblings, key=lambda r: (r.serial or 0))
        index = ordered.index(resonance)
        spellings = nonstereo + stereo
        if index >= len(spellings):
            return None
        atomSetMapping = spellings[index]
    else:
        pool = nonstereo or ambiguous or candidates
        atomSetMapping = sorted(pool, key=_byName)[0]
    name = getattr(atomSetMapping, 'name', None)
    if name:
        return _canonicalAtomName(name, getattr(atomSetMapping, 'mappingType', None))
    return None


def _resonanceIdentity(resonance):
    """Return the NEF identity 4-tuple (chainCode, sequenceCode, residueName, atomName).

    Inverse of the reader's ResonanceGroup naming: an RG named 'A.63' ->
    chain 'A', seq '63'; an RG named without a chain prefix is exported
    under the NEF unassigned-chain convention '@'.  In a native legacy
    project a ResonanceGroup has no name at all - its identity is the
    MolSystem residue it is assigned to, so chain and sequence are then
    taken from ``residue.chain.code`` / ``residue.seqCode`` (the same
    convention as the ``nef_sequence`` rows).  For groups a reader cannot
    place (no name, no residue) or for resonances with no group at all,
    the unassigned-chain form chain '@' + sequence '@<serial>' is written
    (the serial-pinned '@' form the reader's fetchResidueMap resolves to a
    dedicated ResonanceGroup).  Without such a fallback the exported rows
    carry no chain/sequence and the reader drops the assignment on
    reimport.  For a group linked to a residue the standard 3-letter code
    is used (matching the ``nef_sequence`` rows); otherwise the
    importer-assigned group ``ccpCode`` is kept.

    The atom_name is normally ``resonance.name``.  In a native legacy
    project that is None, so the name is recovered from the resonance's
    ``resonanceSet`` atoms matched against the residue's ResidueMapping
    (``_resonanceAtomName``); only resonances with no usable atom identity
    at all are exported in the element@serial pin form."""
    rg = getattr(resonance, 'resonanceGroup', None)
    residue = getattr(rg, 'residue', None) if rg is not None else None
    name = getattr(rg, 'name', None) if rg is not None else None
    if rg is not None and name and '.' in name:
        chainCode, sequenceCode = name.split('.', 1)
    elif rg is not None and name:
        chainCode, sequenceCode = '@', name
    elif residue is not None:
        # Native legacy projects: unnamed groups, identity on the residue
        chain = getattr(residue, 'chain', None)
        chainCode = getattr(chain, 'code', None)
        seqCode = getattr(residue, 'seqCode', None)
        insertCode = (getattr(residue, 'seqInsertCode', None) or '').strip()
        if chainCode is not None and seqCode is not None:
            sequenceCode = (f'{seqCode}{insertCode}') if insertCode else str(seqCode)
        else:
            chainCode, sequenceCode = None, None
    elif rg is not None and getattr(rg, 'serial', None) is not None:
        # Group the reader cannot place: pin it by the group's own serial
        chainCode, sequenceCode = '@', f'@{rg.serial}'
    elif getattr(resonance, 'serial', None) is not None:
        # Resonance without any group: pin it by the resonance serial
        chainCode, sequenceCode = '@', f'@{resonance.serial}'
    else:
        chainCode, sequenceCode = None, None
    if rg is not None:
        residueName = getattr(rg, 'ccpCode', None)
        if residue is not None:
            candidate = _residueName(residue)
            if candidate != 'UNK':
                residueName = candidate
    else:
        residueName = None
    atomName = resonance.name
    if not atomName:
        atomName = _resonanceAtomName(resonance, residue)
    if not atomName:
        # The importer deliberately names resonances None for the reserved
        # 'element@serial' atom form (e.g. H@237); recreate that form so the
        # reimport maps the row to the same resonance
        element, _ = _isotopeToElement(getattr(resonance, 'isotopeCode', None))
        serial = getattr(resonance, 'serial', None)
        if element and serial is not None:
            atomName = f'{element}@{serial}'
    return chainCode, sequenceCode, residueName, atomName


def _isotopeToElement(isotopeCode):
    """Split an isotope code like '13C' into (element, isotopeNumber)."""
    if not isotopeCode or isotopeCode in ('unknown', 'undefined'):
        return None, None
    nn = 0
    while nn < len(isotopeCode) and isotopeCode[nn].isdigit():
        nn += 1
    element = isotopeCode[nn:]
    if not element or not element[0].isalpha():
        return None, None
    return element, (int(isotopeCode[:nn]) if nn else None)


def _objectName(obj):
    """Saveframe suffix name for an object (NEF names must not be '.')."""
    name = getattr(obj, 'name', None)
    if name is None or name in ('', '.'):
        return '1'
    return name


def _uniqueFramecode(db, base):
    """Return a saveframe framecode that is not yet used in db.

    NEF requires unique framecodes per file, but different objects of the
    same class can carry the same name (e.g. two DataSources both named
    'ftt' when spectra are imported from NMRpipe), so a colliding base
    gets a '_2', '_3', ... suffix."""
    framecode = base
    suffix = 2
    while framecode in db:
        framecode = f'{base}_{suffix}'
        suffix += 1
    return framecode


# Key under which the synthesized empty shift-list placeholder frame is
# registered in the per-project framecode map (makeNefDataBlock).
_EMPTY_SHIFT_LIST_KEY = id(None)


def _identityColumns(row, resonances, dimOffset=1, suffixed=True):
    """Fill chain/seq/res/atom columns for a list of (possibly absent) resonances.

    With suffixed=False the unsuffixed single-resonance form is used
    (chemical-shift rows); one resonance only."""
    if not suffixed:
        assert len(resonances) == 1
    for ii, reso in enumerate(resonances):
        suffix = '' if not suffixed else (dimOffset + ii)
        identity = _resonanceIdentity(reso) if reso is not None else (None, None, None, None)
        for tag, value in zip(('chain_code', 'sequence_code', 'residue_name', 'atom_name'),
                             identity):
            row[tag if not suffixed else f'{tag}_{suffix}'] = value


# ---------------------------------------------------------------------------
# saveframe builders
# ---------------------------------------------------------------------------

def _makeMetaData(db):
    sf = db.newSaveFrame('nef_nmr_meta_data', 'nef_nmr_meta_data')
    sf['format_name'] = 'nmr_exchange_format'
    sf['format_version'] = '1.1'
    sf['program_name'] = programName
    sf['program_version'] = programVersion
    timeStamp = commonUtil.getTimeStamp()
    sf['creation_date'] = timeStamp
    sf['uuid'] = commonUtil.getUuid(programName, timeStamp)
    loop = sf.newLoop('nef_program_script', ['program_name', 'script_name'])
    loop.newRow({'program_name': programName, 'script_name': 'exportProject'})
    return sf


def _makeMolecularSystem(db, memopsRoot):
    """nef_molecular_system + nef_sequence from the MolSystem chains."""
    sf = db.newSaveFrame('nef_molecular_system', 'nef_molecular_system')
    columns = ('index', 'chain_code', 'sequence_code', 'residue_name',
               'linking', 'residue_variant', 'cis_peptide')
    loop = sf.newLoop('nef_sequence', columns)

    index = 0
    molSystem = memopsRoot.currentMolSystem
    if molSystem is not None:
        for chain in sorted(molSystem.chains, key=lambda c: (c.code or '')):
            for residue in chain.sortedResidues():
                index += 1
                seqCode = residue.seqCode
                # NB the model default for seqInsertCode is a bare space; strip it
                insertCode = (getattr(residue, 'seqInsertCode', None) or '').strip()
                if seqCode is not None:
                    sequenceCode = (f'{seqCode}{insertCode}') if insertCode else str(seqCode)
                elif insertCode:
                    sequenceCode = insertCode
                else:
                    sequenceCode = '1'
                linking = residue.linking
                if linking == 'none':
                    # NEF v1.1 has no 'none'; an isolated residue is 'single'
                    linking = 'single'
                loop.newRow({
                    'index': index,
                    'chain_code': chain.code,
                    'sequence_code': sequenceCode,
                    'residue_name': _residueName(residue),
                    'linking': linking,
                    'residue_variant': None,
                    'cis_peptide': None,
                })
    return sf


def _shiftLists(nmrProject):
    """Sorted ShiftList objects among the project measurement lists."""
    result = []
    for ml in sorted(nmrProject.measurementLists, key=lambda x: getattr(x, 'serial', 0)):
        if type(ml).__name__ == 'ShiftList':
            result.append(ml)
    return result


def _makeChemicalShiftList(db, shiftList):
    sf = db.newSaveFrame(_uniqueFramecode(db, f'nef_chemical_shift_list_{_objectName(shiftList)}'),
                         'nef_chemical_shift_list')
    sf['atom_chem_shift_units'] = 'ppm'
    if getattr(shiftList, 'details', None):
        sf['ccpn_comment'] = shiftList.details
    sf['ccpn_is_simulated'] = bool(getattr(shiftList, 'isSimulated', False))

    columns = ('chain_code', 'sequence_code', 'residue_name', 'atom_name',
               'value', 'value_uncertainty', 'element', 'isotope_number',
               'ccpn_figure_of_merit', 'ccpn_comment')
    loop = sf.newLoop('nef_chemical_shift', columns)
    # Group shifts into rows.  NB one NEF row may back several Shift objects: an
    # ambiguous atom set (e.g. 'HG%') expands to one resonance per atom set on
    # import, so the shifts of one row carry DIFFERENT name spellings of the
    # same atom.  Writing one row per shift would make such a row expand again
    # on reimport; groups with >1 spelling are therefore written once (in the
    # canonical upper-case '%' form).  Same-spelling duplicates (two identical
    # rows in the source file) are kept as-is, one row per shift.
    groups = {}
    for shift in sorted(shiftList.measurements,
                        key=lambda m: ((getattr(m.resonance, 'serial', 0) or 0),
                                       round(m.value or 0.0, 10))):
        reso = shift.resonance
        element, isotopeNumber = _isotopeToElement(getattr(reso, 'isotopeCode', None))
        row = {
            'value': shift.value,
            'value_uncertainty': shift.error,
            'element': element,
            'isotope_number': isotopeNumber,
        }
        _identityColumns(row, [reso], suffixed=False)
        if getattr(shift, 'figOfMerit', None) is not None:
            row['ccpn_figure_of_merit'] = shift.figOfMerit
        if getattr(shift, 'details', None):
            row['ccpn_comment'] = shift.details
        atom = row['atom_name'] or ''
        key = (row['chain_code'], row['sequence_code'], atom.upper().replace('*', '%'),
               element, isotopeNumber, shift.value)
        group = groups.setdefault(key, {'rows': [], 'spellings': set()})
        group['rows'].append(row)
        group['spellings'].add(atom.upper())
    for group in groups.values():
        rows = group['rows']
        if len(group['spellings']) > 1:
            canonical = [r for r in rows if (r['atom_name'] or '').upper().endswith('%')]
            rows = canonical or rows[:1]
        for row in rows:
            loop.newRow(row)
    return sf


_RESTRRAINT_LIMIT_KEYS = (
    ('target_value', 'targetValue'),
    ('target_value_uncertainty', 'error'),
    ('lower_limit', 'lowerLimit'),
    ('upper_limit', 'upperLimit'),
)


def _constraintLimits(constraint, item=None):
    """Collect the limit/target values from constraint and (dihedral) item."""
    limits = {}
    for tag, attr in _RESTRRAINT_LIMIT_KEYS:
        value = None
        for source in (constraint, item):
            if source is not None and getattr(source, attr, None) is not None:
                value = getattr(source, attr)
                break
        if value is not None:
            limits[tag] = value
    return limits


def _fixedResonance(reso):
    """FixedResonance -> underlying Resonance (for identity export)."""
    underlying = getattr(reso, 'resonance', None)
    return underlying if underlying is not None else reso


def _makeRestraintList(db, constraintList):
    """One NEF saveframe for a constraint list; returns None if unsupported."""
    className = type(constraintList).__name__
    name = _objectName(constraintList)

    if className in ('DistanceConstraintList', 'HBondConstraintList', 'RdcConstraintList',
                     'DihedralConstraintList', 'JCouplingConstraintList',
                     'ChemShiftConstraintList', 'CsaConstraintList'):
        numberResonances = 4 if className.startswith('Dihedral') else (
            1 if className in ('ChemShiftConstraintList', 'CsaConstraintList') else 2)
    else:
        return None

    if className == 'DistanceConstraintList':
        category, loopName = 'nef_distance_restraint_list', 'nef_distance_restraint'
    elif className == 'HBondConstraintList':
        category, loopName = 'nef_distance_restraint_list', 'nef_distance_restraint'
    elif className == 'RdcConstraintList':
        category, loopName = 'nef_rdc_restraint_list', 'nef_rdc_restraint'
    elif className == 'DihedralConstraintList':
        category, loopName = 'nef_dihedral_restraint_list', 'nef_dihedral_restraint'
    else:
        # JCoupling / ChemicalShift / Csa are CCPN-specific restraint types
        category, loopName, restraintType = ('ccpn_restraint_list', 'ccpn_restraint',
                                             className.replace('ConstraintList', '').replace('ChemShift', 'ChemicalShift'))

    sf = db.newSaveFrame(_uniqueFramecode(db, f'{category}_{name}'), category)
    if category == 'ccpn_restraint_list':
        sf['restraint_type'] = restraintType
    else:
        first = list(constraintList.constraints)
        method = getattr(first[0], 'method', None) if first else None
        sf['potential_type'] = method or 'harmonic'
        if className == 'HBondConstraintList':
            sf['restraint_origin'] = 'hbond'

    identityTags = tuple(f'{t}_{ii + 1}' for ii in range(numberResonances)
                         for t in ('chain_code', 'sequence_code', 'residue_name', 'atom_name'))
    columns = ('index', 'restraint_id', 'restraint_combination_id') + identityTags + (
        'weight', 'target_value', 'target_value_uncertainty',
        'lower_linear_limit', 'lower_limit', 'upper_limit', 'upper_linear_limit',
        'ccpn_comment')
    loop = sf.newLoop(loopName, columns)

    index = 0
    for constraint in sorted(constraintList.constraints, key=lambda x: x.serial):
        items = sorted(list(constraint.items), key=lambda x: getattr(x, 'serial', 0))
        if numberResonances == 4:
            # dihedral: the importer creates exactly one item per row, so every
            # item of the constraint gets its own row (resonances on the
            # constraint, limits on the item)
            rowItems = items if items else [None]
        elif numberResonances == 1:
            # single-resonance restraints: resonance on the constraint
            rowItems = [None]
        else:
            # pairwise: one row per constraint; the importer derives the item
            # (product) alternatives from the row's resonances
            rowItems = [items[0]] if items else []

        for item in rowItems:
            index += 1
            if numberResonances == 4:
                resonances = [_fixedResonance(r) for r in constraint.resonances]
                limits = _constraintLimits(constraint, item)
            elif numberResonances == 1:
                reso = getattr(constraint, 'resonance', None)
                resonances = [_fixedResonance(reso)]
                limits = _constraintLimits(constraint)
            else:
                resonances = [_fixedResonance(r) for r in
                              sorted(item.resonances, key=lambda r: (r.serial or 0))]
                limits = _constraintLimits(constraint)

            row = {'index': index, 'restraint_id': constraint.serial}
            _identityColumns(row, resonances)
            row['weight'] = getattr(constraint, 'weight', None)
            row.update(limits)
            if getattr(constraint, 'details', None):
                row['ccpn_comment'] = constraint.details
            loop.newRow(row)
    #
    return sf


def _makeSpectrum(db, nmrProject, experiment, dataSource, shiftListFrameCodes):
    """One nef_nmr_spectrum saveframe for a DataSource (incl. peaks).

    Returns (framecode, peakCodeBySerial) for the peak-restraint links."""

    name = _objectName(dataSource)
    framecode = _uniqueFramecode(db, f'nef_nmr_spectrum_{name}')
    sf = db.newSaveFrame(framecode, 'nef_nmr_spectrum')

    numDim = dataSource.numDim
    sf['num_dimensions'] = numDim

    # Link to a chemical shift list frame (mandatory field).  When the
    # project has no shift lists at all, makeNefDataBlock registers one
    # synthesized empty placeholder frame (the importer reads it as an
    # empty ShiftList, which is harmless), shared by every spectrum.
    expShiftList = getattr(experiment, 'shiftList', None)
    if expShiftList is not None and id(expShiftList) in shiftListFrameCodes:
        sf['chemical_shift_list'] = shiftListFrameCodes[id(expShiftList)]
    else:
        sf['chemical_shift_list'] = next(iter(shiftListFrameCodes.values()))

    if experiment.name:
        sf['experiment_type'] = experiment.name
    refExperiment = getattr(experiment, 'refExperiment', None)
    if refExperiment is not None and getattr(refExperiment, 'name', None):
        sf['experiment_classification'] = refExperiment.name

    # A DataSource linked to a data file carries the link in the NEF so a
    # plain reimport (Load NEF, no relink step) restores it: the importer
    # (v2io.NefIo) reads ccpn_spectrum_file_path + the ccpn_file_* items
    # into a dataStore (addDataStore) and point_count/total_point_count
    # onto the data dimensions.  Unlinked DataSources export exactly as
    # before (same columns, no file items).
    dataStore = getattr(dataSource, 'dataStore', None)

    # ------------------------------------------------------------------ dims
    dimColumns = ('dimension_id', 'axis_unit', 'axis_code', 'spectrometer_frequency',
                  'spectral_width', 'value_first_point', 'folding', 'is_acquisition')
    dimLoop = sf.newLoop('nef_spectrum_dimension', dimColumns)
    for dataDim in sorted(dataSource.dataDims, key=lambda x: x.dim):
        dim = dataDim.dim
        expDim = experiment.findFirstExpDim(dim=dim)
        expDimRef = expDim.findFirstExpDimRef() if expDim is not None else None
        dataDimRef = dataDim.findFirstDataDimRef()

        axisCode = None
        if expDimRef is not None and getattr(expDimRef, 'isotopeCodes', None):
            axisCode = expDimRef.isotopeCodes[0]
        row = {'dimension_id': dim,
               'axis_unit': (expDimRef.unit if expDimRef is not None else None) or 'ppm',
               'axis_code': axisCode or '1H'}
        if expDimRef is not None:
            if expDimRef.sf is not None:
                row['spectrometer_frequency'] = expDimRef.sf
            if expDimRef.isFolded:
                row['folding'] = 'mirror'
            else:
                row['folding'] = 'circular'
            if (expDimRef.sf and dataDim.numPoints and dataDim.valuePerPoint):
                row['spectral_width'] = dataDim.valuePerPoint * dataDim.numPoints / expDimRef.sf
        if dataDimRef is not None and dataDimRef.refValue is not None:
            row['value_first_point'] = dataDimRef.refValue
        row['is_acquisition'] = bool(expDim.isAcquisition) if expDim is not None else False
        dimLoop.newRow(row)

    # ------------------------------------------------------------- file link
    if dataStore is not None:
        # absolute path: the importer hands it straight to addDataStore
        sf['ccpn_spectrum_file_path'] = dataStore.fullPath
        for nefKey, attr in (
            ('ccpn_file_type', 'fileType'),
            ('ccpn_file_header_size', 'headerSize'),
            ('ccpn_file_byte_number', 'nByte'),
            ('ccpn_file_number_type', 'numberType'),
            ('ccpn_file_is_big_endian', 'isBigEndian'),
            ('ccpn_file_complex_stored_by', 'complexStoredBy'),
        ):
            value = getattr(dataStore, attr, None)
            if value is not None:
                sf[nefKey] = value
        # the point counts (and the matrix block layout) ride on the ccpn
        # extension loop - the importer reads them from
        # ccpn_spectrum_dimension (NOT from nef_spectrum_dimension, whose
        # column set it does not extend).  dimension_block_size is load
        # bearing: without it the importer guesses block sizes from the
        # grid (determineBlockSizes) and misreads NMRpipe files whose real
        # layout is (npts_dim1, 1, ...) - e.g. (427, 1) read as (128, 32).
        pcLoop = sf.newLoop('ccpn_spectrum_dimension',
                            ('dimension_id', 'point_count', 'total_point_count',
                             'dimension_block_size'))
        storeBlocks = tuple(getattr(dataStore, 'blockSizes', None) or ())
        for dataDim in sorted(dataSource.dataDims, key=lambda x: x.dim):
            row = {'dimension_id': dataDim.dim,
                   'point_count': dataDim.numPoints,
                   'total_point_count': dataDim.numPointsOrig}
            ii = dataDim.dim - 1
            if 0 <= ii < len(storeBlocks):
                row['dimension_block_size'] = int(storeBlocks[ii])
            pcLoop.newRow(row)

    # ---------------------------------------------------------------- transfer
    transferColumnSet = ('dimension_1', 'dimension_2', 'transfer_type', 'is_indirect')
    transferLoop = sf.newLoop('nef_spectrum_dimension_transfer', transferColumnSet)
    for transfer in sorted(experiment.expTransfers,
                           key=lambda x: tuple(sorted(r.expDim.dim for r in x.expDimRefs))):
        dims = sorted(expDimRef.expDim.dim for expDimRef in transfer.expDimRefs)
        transferLoop.newRow({'dimension_1': dims[0], 'dimension_2': dims[1],
                             'transfer_type': transfer.transferType,
                             'is_indirect': not transfer.isDirect})

    # ---------------------------------------------------------------- peaklist
    peakLists = [pl for pl in sorted(dataSource.peakLists,
                                     key=lambda x: getattr(x, 'serial', 0))]
    peakList = peakLists[0] if peakLists else None
    if peakList is not None:
        if peakList.name:
            sf['ccpn_peaklist_name'] = peakList.name
        if getattr(peakList, 'details', None):
            sf['ccpn_peaklist_comment'] = peakList.details
        sf['ccpn_peaklist_is_simulated'] = bool(getattr(peakList, 'isSimulated', False))

    # ---------------------------------------------------------------- peaks
    dimCount = numDim
    peakColumns = tuple(['index', 'peak_id', 'volume', 'volume_uncertainty',
                         'height', 'height_uncertainty'])
    for dim in range(1, dimCount + 1):
        peakColumns += (f'position_{dim}', f'position_uncertainty_{dim}')
    for dim in range(1, dimCount + 1):
        peakColumns += (f'chain_code_{dim}', f'sequence_code_{dim}',
                        f'residue_name_{dim}', f'atom_name_{dim}')
    peakColumns += ('ccpn_annotation', 'ccpn_comment', 'ccpn_figure_of_merit')
    peakLoop = sf.newLoop('nef_peak', peakColumns)

    peakCodeBySerial = {}
    if peakList is not None:
        for ii, peak in enumerate(sorted(peakList.peaks, key=lambda x: x.serial)):
            peakCodeBySerial[peak.serial] = framecode
            base = {'index': ii + 1, 'peak_id': peak.serial,
                    'volume': peak.volume, 'volume_uncertainty': None,
                    'height': peak.height, 'height_uncertainty': None}
            for pd in peak.sortedPeakDims():
                base[f'position_{pd.dim}'] = pd.value
                base[f'position_uncertainty_{pd.dim}'] = pd.valueError
            attribs = {}
            if getattr(peak, 'annotation', None):
                attribs['ccpn_annotation'] = peak.annotation
            if getattr(peak, 'details', None):
                attribs['ccpn_comment'] = peak.details
            if getattr(peak, 'figOfMerit', None) is not None:
                attribs['ccpn_figure_of_merit'] = peak.figOfMerit

            contribs = list(peak.peakContribs)
            if not contribs:
                row = dict(base)
                _identityColumns(row, [None] * dimCount)
                row.update(attribs)
                peakLoop.newRow(row)
                continue
            # One row per alternative assignment (the importer's multi-row form)
            for contrib in contribs:
                row = dict(base)
                contribDimMap = {}
                for pdc in contrib.peakDimContribs:
                    contribDimMap[pdc.peakDim.dim] = pdc.resonance
                _identityColumns(row, [contribDimMap.get(dim) for dim in range(1, dimCount + 1)])
                row.update(attribs)
                peakLoop.newRow(row)
    #
    return framecode, peakCodeBySerial


class _EmptyShiftList:
    """Placeholder for the mandatory shift-list cross-reference on projects
    that have no ShiftLists."""
    name = None
    details = None
    isSimulated = False
    measurements = []


def _makePeakRestraintLinks(db, links):
    sf = db.newSaveFrame('nef_peak_restraint_links', 'nef_peak_restraint_links')
    columns = ('nmr_spectrum_id', 'peak_id', 'restraint_list_id', 'restraint_id')
    loop = sf.newLoop('nef_peak_restraint_link', columns)
    for spectrumCode, peakSerial, restraintCode, restraintSerial in links:
        loop.newRow({'nmr_spectrum_id': spectrumCode, 'peak_id': peakSerial,
                     'restraint_list_id': restraintCode, 'restraint_id': restraintSerial})
    return sf


# ---------------------------------------------------------------------------
# top level
# ---------------------------------------------------------------------------

def makeNefDataBlock(memopsRoot):
    """Build a StarIo.NmrDataBlock (contemporary NEF v1.1) from a legacy project."""
    nmrProject = memopsRoot.currentNmrProject
    if nmrProject is None:
        raise ValueError(f'no NmrProject found under memopsRoot {memopsRoot}')

    db = StarIo.NmrDataBlock()
    db.name = getattr(memopsRoot, 'name', None) or 'nmr_project'

    _makeMetaData(db)
    _makeMolecularSystem(db, memopsRoot)

    shiftListFrameCodes = {}
    for shiftList in _shiftLists(nmrProject):
        sf = _makeChemicalShiftList(db, shiftList)
        shiftListFrameCodes[id(shiftList)] = sf['sf_framecode']
    if not shiftListFrameCodes:
        # Project without shift lists: synthesize one empty list frame for
        # the mandatory cross-reference on the spectrum saveframes (the
        # importer reads it as an empty ShiftList, which is harmless).
        # Created once and shared by all - framecodes must be unique
        # per file.
        sf = _makeChemicalShiftList(db, _EmptyShiftList())
        shiftListFrameCodes[_EMPTY_SHIFT_LIST_KEY] = sf['sf_framecode']

    restraintFrameCodes = {}
    for store in sorted(nmrProject.nmrConstraintStores, key=lambda x: getattr(x, 'serial', 0)):
        for constraintList in sorted(store.constraintLists, key=lambda x: getattr(x, 'serial', 0)):
            sf = _makeRestraintList(db, constraintList)
            if sf is not None:
                restraintFrameCodes[id(constraintList)] = sf['sf_framecode']

    spectrumFrameByDataSource = {}
    for experiment in sorted(nmrProject.experiments, key=lambda x: getattr(x, 'serial', 0)):
        for dataSource in sorted(experiment.dataSources, key=lambda x: getattr(x, 'serial', 0)):
            framecode, _ = _makeSpectrum(db, nmrProject, experiment, dataSource,
                                         shiftListFrameCodes)
            spectrumFrameByDataSource[id(dataSource)] = framecode

    # peak <--> restraint links
    links = []
    for store in sorted(nmrProject.nmrConstraintStores, key=lambda x: getattr(x, 'serial', 0)):
        for constraintList in sorted(store.constraintLists, key=lambda x: getattr(x, 'serial', 0)):
            restraintCode = restraintFrameCodes.get(id(constraintList))
            if restraintCode is None:
                continue
            for constraint in constraintList.constraints:
                peaks = getattr(constraint, 'peaks', None)
                if not peaks:
                    continue
                for peak in peaks:
                    peakList = peak.peakList
                    dataSource = peakList.dataSource if peakList is not None else None
                    spectrumCode = spectrumFrameByDataSource.get(id(dataSource)) \
                        if dataSource is not None else None
                    if spectrumCode is None:
                        continue
                    links.append((spectrumCode, peak.serial, restraintCode, constraint.serial))
    if links:
        _makePeakRestraintLinks(db, links)

    return db


def exportProject(memopsRoot, fileName):
    """Write the NEF (contemporary, v1.1) file for memopsRoot to fileName.

    Returns fileName."""
    dataBlock = makeNefDataBlock(memopsRoot)
    text = dataBlock.toString()
    dirName = os.path.dirname(os.path.abspath(fileName))
    if not os.path.isdir(dirName):
        os.makedirs(dirName)
    with open(fileName, 'w') as fp:
        fp.write(text)
    return fileName


if __name__ == '__main__':
    # quick CLI: exportProject <memopsRoot-directory-or-project-dir> <out.nef>
    # (the full command-line interface is the ccpnmr-nef console entry
    # point in ccpnmr/nefCli.py)
    if len(sys.argv) != 3:
        sys.stderr.write(f'usage: {sys.argv[0]} <project-directory> <output.nef>\n')
        sys.exit(1)
    from memops.general import Io as memopsIo
    root = memopsIo.loadProject(sys.argv[1])
    exportProject(root, sys.argv[2])
    print('wrote', sys.argv[2])
