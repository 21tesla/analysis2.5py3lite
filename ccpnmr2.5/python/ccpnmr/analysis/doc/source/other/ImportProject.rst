================
Project : Import
================

**Import information from a non-CCPN source into the current project**

This submenu allows the user to add data to a CCPN project using data stored in
file formats external to the CCPN system. Typically this is used to add
ancilliary and derived NMR data such as chemical shift lists, peak lists and
macromolecular coordinates etc. The data matrices for spectra sould be imported
via the `Open Spectra`_ option, *not* by using these functions.

**Direct Import Options**

The user may select an option to load a directly supported format, after which
the user will be presented with a file browser to locate the file to be loaded. 

BioMagResBank data may be imported by selecting the option for the relevant version of
NMR-STAR; the textual NMR depostion format. Virtually everything present in the
NMR-STAR file will be imported, i.e. any chemical shift lists, molecular information,
peak lists, citations, sample preparation etc. Because the BMRB data is normally
considered as part of an NMR database deposition the imported information will be
linked into an "Entry" record that may be viewed via the ECI_ option.

Well-formatted, official PDB version 3.20 formated files containing structural
coordinate data and header information may be loaded into the CCPN project by
seleting the "PDB 3.20" option. Using this option not only enters atomic
coordinate data into the CCPN project but also all of the header information that
was specified when the PDB entry was deposited; sample specification, citations etc.
Such deposition data may be viewed via the ECI_ option.

Non-standard PDB-style files may be loaded with the "Coordinates" option. This does
not consider any header information, only the ATOM records, and thus only loads the
main atomic coordinates. The advantage to using this option is that it is fairly
flexible as to which convention was used for the names of atoms, and this need not be
known in advance. The loader attempts to automatically deduce the best naming
convention and uses fall-back atom names if anything still does not fit as expected.

**Importing from other formats**

The Format Converter option should be used for or importing data from any formats
not directly listed in the menu. Please note that Format Converter documentation
is not currently covered in these pages.

.. _`Open Spectra`: ../popups/OpenSpectrumPopup.html

.. _ECI: ../popups/EntryCompletionPopup.html
