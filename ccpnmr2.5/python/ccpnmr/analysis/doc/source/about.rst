.. CcpNmr Python API Programming Guide

=====================
About CcpNmr Analysis
=====================

CcpNmr Analysis is a graphics-based interactive NMR spectrum visualisation, resonance assigment
and data analysis program. CcpNmr Analysis can be considered a platform for almost all the NMR
data described by the CCPN `data model`_ and a place from which to interact with connected
non-CCPN programs, for example those integrated in the Extend-NMR_ project. To name only a subset
of its capabilities, you can use CcpNmr Analysis to assign atomic resonances (proteins, nucleic
acid, glycans & small molecules), measure relaxation rates, measure kinetic paremeters, calculate
distance constraints, co-ordinate structure generation, edit molecular information and as a
starting point to develop new software and algorithms for NMR.

**Development**

CcpNmr Analysis was primarily witten by Wayne Boucher and Tim Stevens at the University of
Cambridge. Most of the computer code for Analysis is written in Python. Only speed-critical
functions are performed by code written in the fast, compiled language C. Such C functions include
the calculation of contours and mathematically intensive algorithms. The Python part of the
program consists of a series of integrated graphical windows ("Popups") and an underlying layer of
`Python library functions`_. The graphical elements allow the user to enter information and to
view the status of the data, while the library functions manipulate the CCPN `data model`_ objects
to record the scientific information.

CcpNmr Analysis has its own CCPN `data model`_ package, which means it uses CCPN technology to create a
program-specific part of the Python API and thus allows program information (colours, window
positions etc.) to be recorded and stored as XML files. 

**Objectives**

At the start of the CCPN project the requirement was for an NMR data analysis program that used a
modern graphical user interface and could run on many types of computer. It would be supported and
maintained by CCPN and would allow modification and extension, including for new NMR techniques.
The first version of Analysis was released in 2005. Analysis is built directly on the CCPN data
model and its design is partly inspired by the older ANSIG[1]. and SPARKY[2] programs, but it has
continued to develop from the suggestions, requirements and computational contributions of its
user community. Analysis is freely available to academic and non-profit institutions. Commercial
users are required to subscribe to CCPN for a moderate fee.

Main Program Concepts
=====================

If you've not seen this documentation before you may want to know what it's all about, so here are
some highlights of the CcpNmr Analysis software. 

Spectrum Windows
----------------

An Analysis project may contain an almost limitless number of spectrum windows. The windows are
are inherently N-dimensional with scrollbars for not only the screen dimensions but also for
orthogonal planes, with the ability to select any plane thickness. A window can be divided into
several strips for easy comparison of different regions of spectra. Many spectra may be
superimposed in the same window where their contours and peaks are readily toggled on or off.
Navigation is achieved by using the mouse or keyboard and there are inbuilt navigation functions
to easily find orthogonal planes and return-peak positions etc. Many functions may be applied to
crosspeaks directly from the window menu. For example peaks may be assigned, deleted, unaliased
and shift matched. Several of these functions can be used on several peaks at once to improve user
efficiency. For example, columns of NOE peaks derived from the same amide resonances may be
assigned to this amide at the same time.

Molecules
---------

Polymer chains and small molecules are readily put into NMR projects. Sequences may be imported
from file or entered directly by the user and from this Analysis will build the molecules with all
of their NMR assignable atoms. Many molecules of different types can be included and may be
connected together into chains. For example, a GIP-anchored glycoprotein may be constructed by
joining protein, sugar and lipid components. By using data provided by the PDBe at the EBI, CcpNmr
software has access to a large number small molecule templates - those that have appeared in PDB
structures.

Tables
------

Virtually all of the information within a project is available to the user via a graphical
interface (and a Python shell should you be brave enough to use it) and much of the commonly used
information is presented in tabular form. These tables are used to display peak lists, chemical
shifts, constraints, coordinates, spectrum configuration and the like. To allow the user to change
information (peak position, contour colour, experiment name to name only a few...) they often have
editable columns. The rows of the table may be sorted on any of the column types, may be filtered
according to a search expression and may be selected (often several at once) to apply specific
functions. Also, the data in a table may be exported to a text file, output as PostScript and if
numerical may be plotted in a graph.

Resonance Assignment
--------------------

Assignment in CcpNmr software is a two-step process proceeding via an intermediate Resonance
object. This allows the user to represent anonymous but connected assignment states, and allows
atomic assignment to be made to several peaks at once. Most crosspeak assignment is made by the
user choosing a resonance (which need not be assigned to specific atoms) from a curated/ranked
list. The choice is made with a single click (and is readily reversible) from a list of
possibilities that are close in chemical shift. Structural information can also be used in
assignment. Here through-space linked resonances may be ranked according to their distance in an
intermediate structure.

Assignment of resonances to specific atoms is achieved by selecting the atom on a display showing
the well-curated molecular information. This needs to be done only once for each resonance as
all peaks which correspond to the same resonance will automatically share the atom information.
A resonance may be assigned, where appropriate, in a stereospecific or non-stereospecific manner.
For example, it is possible to say that two peaks represent two different hydrogen beta atoms in a
residue, with different chemical shifts, but without necessarily specifying the stereochemical
arrangement of the two atoms.

There are many tools designed to expedite the resonance assignment process, including using root
resonance locations (e.g. amide) to direct peak picking and assignment in higher dimensionality
spectra, automated matching of peak positions for sequential protein backbone assignment and
chemical shift plus structure based filtering of NOE/through-space assignment possibilities.

Structure Generation
--------------------

Analysis can be used to generate distance and dihedral angle restraints for structure generation.
Distance restraints may be generated from assigned NOESY peaks, or may be created by performing
shift matching on unassigned peaks. The potentially ambiguous constraints thus output may then be
used by programs, such as ARIA, which are able to take input data from a CCPN project and write the
results back. The violations that result from a structure generation cycle may be imported into Analysis, from
where the user can readily follow a link to the peaks which were used in the generation of the
violated constraint. Connections to the CING software allow easy validation of macromolecular
structures.

Data Analyses
-------------

Analysis has specialist tools to extract relaxation rates, chemical shift changes, kinetic
parameters, scalar couplings etc. These are designed to make such tasks in NMR less tedious, and
the program aims to automate as much of the simple parts of the processes as possible. For example
when following chemical shift changes in titration experiments Analysis automatically tracks the
trajectories of shifting peaks so that the titration points can be considered as an analytic
group. The program also goes on to fit equation curves to the data and extrct the relevant
paremeters.

Reference Information
---------------------

All the CcpNmr programs have access to a library of reference information. This includes chemical
compound descriptions, chemical shift distributions, isotope information, idealised residue
coordinates etc. This is often used implicitly within Analysis, so that the user doesn't have to
worry about how to get hold of such information. Some of the data is visualised where it can be
helpful. For example, chemical shift distributions during assignment.

.. rubric:: Footnotes

.. [1] P.J. Kraulis, "ANSIG: A Program for the Assignment of Protein 1H 2D NMR
       spectra by Interactive Graphics" (1989) J. Magn. Reson 24, pp 627-633

.. [2] T.D. Goddard and D.G. Kneller, SPARKY 3, University of California,
       San Francisco

.. _`data Model`: http://www.ccpn.ac.uk/api-documentation/ccpnmr/ccpnmr2.0/python/doc/api.html
.. _Extend-NMR: http://www.extend-nmr.eu/
.. _`Python library functions`: codingLibrary.html
