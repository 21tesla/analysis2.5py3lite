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

def _format_tab_atom(atom):
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
    chain = residue.chain
    chain_code = chain.code if chain else None
    
    if chain_code:
        return f"{chain_code}-{one_letter}{seq_num}-{atom.name}"
    else:
        return f"{one_letter}{seq_num}-{atom.name}"

def get_tab_ass(peak):
    # Try to find a proton atom first
    for pd in peak.sortedPeakDims():
        atom = _getPeakDimAtom(pd)
        if atom is not None and atom.chemAtom is not None and atom.chemAtom.elementSymbol == 'H':
            return _format_tab_atom(atom)
    # Fallback to any atom
    for pd in peak.sortedPeakDims():
        atom = _getPeakDimAtom(pd)
        if atom is not None:
            return _format_tab_atom(atom)
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
