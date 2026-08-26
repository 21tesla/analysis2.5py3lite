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
  * DataStore / raw spectrum matrix data is not part of NEF and is not
    exported (NEF carries peaks + shifts + restraints + metadata).

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


def _resonanceIdentity(resonance):
    """Return the NEF identity 4-tuple (chainCode, sequenceCode, residueName, atomName).

    Inverse of the reader's ResonanceGroup naming: RG name 'A.63' ->
    chain 'A', seq '63'; an RG named without a chain prefix is exported
    under the NEF unassigned-chain convention '@'.  For a group linked to a
    MolResidue the standard 3-letter code is used (matching the
    ``nef_sequence`` rows); otherwise the importer-assigned group
    ``ccpCode`` is kept."""
    rg = getattr(resonance, 'resonanceGroup', None)
    if rg is not None:
        name = rg.name
        if name and '.' in name:
            chainCode, sequenceCode = name.split('.', 1)
        elif name:
            chainCode, sequenceCode = '@', name
        else:
            chainCode, sequenceCode = None, None
        residueName = getattr(rg, 'ccpCode', None)
        residue = getattr(rg, 'residue', None)
        if residue is not None:
            candidate = _residueName(residue)
            if candidate != 'UNK':
                residueName = candidate
    else:
        chainCode = sequenceCode = residueName = None
    atomName = resonance.name
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
        for chain in sorted(molSystem.chains):
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
    sf = db.newSaveFrame(f'nef_chemical_shift_list_{_objectName(shiftList)}',
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
    for shift in sorted(shiftList.measurements):
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

    sf = db.newSaveFrame(f'{category}_{name}', category)
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
                resonances = [_fixedResonance(r) for r in sorted(item.resonances)]
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
    framecode = f'nef_nmr_spectrum_{name}'
    sf = db.newSaveFrame(framecode, 'nef_nmr_spectrum')

    numDim = dataSource.numDim
    sf['num_dimensions'] = numDim

    # Link to a chemical shift list frame (mandatory field)
    expShiftList = getattr(experiment, 'shiftList', None)
    if expShiftList is not None and id(expShiftList) in shiftListFrameCodes:
        sf['chemical_shift_list'] = shiftListFrameCodes[id(expShiftList)]
    elif shiftListFrameCodes:
        sf['chemical_shift_list'] = next(iter(shiftListFrameCodes.values()))
    else:
        # Project without shift lists: synthesize an empty list frame for
        # the mandatory cross-reference (importer will read it as an empty
        # ShiftList, which is harmless).
        empty = _makeChemicalShiftList(db, _EmptyShiftList())
        sf['chemical_shift_list'] = empty['sf_framecode']

    if experiment.name:
        sf['experiment_type'] = experiment.name
    refExperiment = getattr(experiment, 'refExperiment', None)
    if refExperiment is not None and getattr(refExperiment, 'name', None):
        sf['experiment_classification'] = refExperiment.name

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
    if len(sys.argv) != 3:
        sys.stderr.write(f'usage: {sys.argv[0]} <project-directory> <output.nef>\n')
        sys.exit(1)
    from memops.general import Io as memopsIo
    memopsIo.loadProject(sys.argv[1])
    root = memopsIo.memopsRoot
    exportProject(root, sys.argv[2])
    print('wrote', sys.argv[2])
