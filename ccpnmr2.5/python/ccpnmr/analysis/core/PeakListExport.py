"""Core module for exporting peak lists to external formats (.tab, .nef, .peaks).
"""

import os
from ccpnmr.analysis.core import PeakBasic

def _getPeakDimAtom(peak_dim):
    for contrib in peak_dim.peakDimContribs:
        resonance = contrib.resonance
        resonanceSet = getattr(resonance, 'resonanceSet', None)
        if resonanceSet is not None:
            atomSet = resonanceSet.findFirstAtomSet()
            if atomSet is not None:
                atom = atomSet.findFirstAtom()
                if atom is not None:
                    return atom
    return None

def _format_tab_atom_new(atom):
    residue = atom.residue
    res_code = residue.ccpCode or ""
    CCPN_TO_ONE_LETTER = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    one_letter = CCPN_TO_ONE_LETTER.get(res_code, res_code[0] if res_code else "?")
    seq_num = residue.seqCode
    return f"{one_letter}{seq_num}{atom.name}"

def get_tab_ass(peak):
    parts = []
    for pd in peak.sortedPeakDims():
        atom = _getPeakDimAtom(pd)
        if atom is not None:
            parts.append(_format_tab_atom_new(atom))
        else:
            parts.append(None)
            
    if all(p is not None for p in parts):
        return "-".join(parts)
        
    # If not all dimensions are assigned, fallback to any assigned atom using old style format, or "*"
    for pd in peak.sortedPeakDims():
        atom = _getPeakDimAtom(pd)
        if atom is not None:
            res_code = atom.residue.ccpCode or ""
            CCPN_TO_ONE_LETTER = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            one_letter = CCPN_TO_ONE_LETTER.get(res_code, res_code[0] if res_code else "?")
            return f"{one_letter}{atom.residue.seqCode}-{atom.name}"
    return "*"

def exportTabPeaks(peakList, filePath):
    """
    Export the peaks of a CCPN PeakList into an NMRdraw .tab file.
    """
    dataSource = peakList.dataSource
    num_dim = dataSource.numDim
    axis_letters = "XYZWUV"[:num_dim]
    
    # Header lines
    header_vars = ["INDEX"]
    for letter in axis_letters:
        header_vars.append(f"{letter}_PPM")
    header_vars += ["HEIGHT", "VOL", "ASS"]
    
    vars_line = "VARS   " + " ".join(header_vars)
    
    # Format line
    format_parts = ["%5d"]
    for _ in range(num_dim):
        format_parts.append("%8.3f")
    format_parts += ["%+e", "%+e", "%s"]
    format_line = "FORMAT " + " ".join(format_parts)
    
    lines = [
        vars_line,
        format_line,
        "",
        "NULLVALUE -666",
        "NULLSTRING *",
        ""
    ]
    
    for ii, peak in enumerate(sorted(peakList.peaks, key=lambda x: x.serial)):
        row_id = ii + 1
        pos_parts = []
        for pd in peak.sortedPeakDims():
            pos_parts.append(pd.value if pd.value is not None else 0.0)
            
        height = peak.height if peak.height is not None else 0.0
        volume = peak.volume if peak.volume is not None else 0.0
        ass = get_tab_ass(peak)
        
        # Build row tokens
        tokens = [f"{row_id:5d}"]
        for pos in pos_parts:
            tokens.append(f"{pos:8.3f}")
        tokens.append(f"{height:+e}")
        tokens.append(f"{volume:+e}")
        tokens.append(ass)
        
        lines.append(" ".join(tokens))
        
    with open(filePath, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def exportNefPeaks(peakList, filePath):
    """
    Export the peaks of a CCPN PeakList into a NEF (.nef) file.
    """
    from ccpnmr import nefExport
    nefExport.exportProject(peakList.root, filePath)


def _isotope_to_xeasy_axis(isotope_code):
    if not isotope_code:
        return "H"
    isotope_upper = isotope_code.upper()
    if "H" in isotope_upper:
        return "H"
    if "N" in isotope_upper:
        return "N"
    if "C" in isotope_upper:
        return "C"
    return "H"


def _get_xeasy_dim_ass(peak_dim):
    atom = _getPeakDimAtom(peak_dim)
    if atom is not None:
        return f"{atom.name}.{atom.residue.seqCode}"
    return "-"


def exportXeasyPeaks(peakList, filePath):
    """
    Export the peaks of a CCPN PeakList into an XEASY .peaks file.
    """
    dataSource = peakList.dataSource
    num_dim = dataSource.numDim
    dataDims = sorted(dataSource.dataDims, key=lambda x: x.dim)
    
    # Get isotope codes
    dimIsotopes = {}
    for dataDim in dataDims:
        ref = dataDim.findFirstDataDimRef()
        if ref is not None and ref.expDimRef is not None and ref.expDimRef.isotopeCodes:
            dimIsotopes[dataDim.dim] = ref.expDimRef.isotopeCodes[0]
        else:
            dimIsotopes[dataDim.dim] = "unknown"
            
    xeasy_axes = [_isotope_to_xeasy_axis(dimIsotopes[dd.dim]) for dd in dataDims]
    
    lines = [
        f"# Number of dimensions {num_dim}",
        f"#FORMAT xeasy{num_dim}D"
    ]
    for dim_idx, axis in enumerate(xeasy_axes):
        lines.append(f"#INAME {dim_idx + 1} {axis}")
        
    spec_name = dataSource.name or "Spectrum"
    lines.append(f"#SPECTRUM {spec_name} " + " ".join(xeasy_axes))
    
    for ii, peak in enumerate(sorted(peakList.peaks, key=lambda x: x.serial)):
        row_id = ii + 1
        pos_parts = []
        assignments = []
        for pd in peak.sortedPeakDims():
            pos_parts.append(pd.value if pd.value is not None else 0.0)
            assignments.append(_get_xeasy_dim_ass(pd))
            
        # Get height or volume
        volume = peak.volume if peak.volume is not None else 0.0
        volume_err = 0.0
        
        # Build line
        line = f"{row_id:7d}"
        for pos in pos_parts:
            line += f" {pos:8.3f}"
        line += f" 1 U {volume:11.3E} {volume_err:9.3E} e 0"
        for ass in assignments:
            line += f" {ass:<9s}"
        lines.append(line)
        
    with open(filePath, "w") as fp:
        fp.write("\n".join(lines) + "\n")
