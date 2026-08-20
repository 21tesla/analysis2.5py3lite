
"""
======================COPYRIGHT/LICENSE START==========================

EntryCompletionFrame.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.

A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""

import os
import time
import tarfile

from memops.api import Implementation

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showYesNo, showOkCancel, showWarning, showError
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.Text import Text

from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect import FileType

from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport
from ccpnmr.format.converters.PdbFormat import PdbFormat
from ccpnmr.format.converters.CnsFormat import CnsFormat

from ccpnmr.eci.EciShiftAnalysis import analyseChemicalShifts

from memops.general.Io import getUserDataPath

#from ccpnmr.analysis.core.QualityControlBasic import analyseChemicalShifts

# # # # #  F O R   L A T E R  # # # # #
# More automated attribs for Nmr list details
# Shift notifiers for PPM check - or recalculate button
# Derive more global data ???

# Missing data rows can be clicked to take you to the right slot

# Buttons that will tar-gzip your project and send them to AutoDep
#   for projects with structures and BMRB for those just with NMR data.

# TODO/BUGS for current software:
#   Citation - author pulldown fonts change

# API name : Human Name : Documentation
NMR_LISTS_DATA = [('ShiftList','Shift lists','Assigned NMR chemical shifts'),
                  ('PeakList','Peak lists','Spectral peak lists'),
                  ('JCouplingList','J Coupling Lists','Lists of J couplings'),
                  ('ShiftDifferenceList','Shift difference lists','Chemical shift difference lists'),
                  ('ShiftAnisotropyList','CSA lists','Chemical Shift anisotropy lists'),
                  ('NoeList','Hetero/Homo NOE lists','Heteronuclear and homonuclear NOE lists'),
                  ('T1List','T1 lists','T1 (R1) NMR relaxation data'),
                  ('T2List','T2 lists','T2 (R2) NMR relaxation data'),
                  ('T1rhoList','T1rho lists','T1rho (R1rho) NMR relaxation data'),
                  ('RdcList','RDC lists','Residual dipolar coupling lists'),
                  ('HExchRateList','H exchange lists','Hydrogen exchange rates data'),
                  ('HExchProtectionList','H protection lists','Hydrogen exchange protection data'),
                  ('DipolarRelaxList','Dipole-dipole relax lists','Dipole-dipole relaxation data'),
                  ('IsotropicS2List','Order params','Isotropic S^2 Order parameters'),
                  ('SpectralDensityList','Spectral density','Spectral density factors'),
                  ('PkaList','pKa lists','NMR-derived pKas (acid dissociation constants)'),
                  ('DataList','Other data lists','Other NMR-derived data lists'),]

#['Dipole-dipole cross corrs','Dipole-dipole cross correlation data']
#['Dipole-CSA cross corrs','Dipole-CSA cross correlation data']
#['pH transitions','NMR-derived pH transitions (pKa, pHmid)']
#['D/H fractionation','NMR-derived D/H fractionation factors']
#['Calc shift lists','Theoretical (calculated) chemical shift values']

STANDARD_ISOTOPES = ['1H','13C','15N','31P','2H','29Si','19F','17O']

# Lots of things below could be derived from the data model!
MEASUREMENT_LIST_CLASSES = ['DipolarRelaxList', 'HExchProtectionList',
                            'HExchRateList', 'JCouplingList', 'NoeList',
                            'RdcList', 'ShiftAnisotropyList',
                            'ShiftDifferenceList', 'ShiftList',
                            'T1List', 'T1rhoList', 'T2List']

DERIVED_LIST_CLASSES =  ['DataList', 'IsotropicS2List',
                         'PkaList', 'SpectralDensityList']

CONSTRAINT_NAME_DATA = {'ChemShiftConstraintList': 'Chemical Shift Constraints',
                        'CsaConstraintList':       'CSA Constraints',
                        'DihedralConstraintList':  'Dihedral Angle Constraints',
                        'DistanceConstraintList':  'Distance (NOE) Constraints',
                        'HBondConstraintList':     'Hydrogen Bond Constraints',
                        'JCouplingConstraintList': 'Coupling Constant Constraints',
                        'RdcConstraintList':       'RDC Constraints'}

BOOL = ['yes', 'no']

UNDEFINED = 'undefined'

TITLES = ['Dr', 'Prof', 'Mr', 'Ms', 'Mrs', 'Miss', 'Sir', 'Dame'] 

FAM_TITLES = ['Jr.', 'Sr.', 'III', 'IV', 'V'] 

POS_ROLES = ['principal investigator', 'responsible scientist',
             'investigator']

ORGN_TYPES = ['academic', 'commercial', 'government', 'other']

VOLUME_UNITS = ['ul','ml','l']

# liquid and solid not in data model, but is in NMR-STAR

EXP_STATES  = ['liquid', 'solid', 'isotropic', 'anisotropic', 'ordered', 'powder', 'crystal']

SAMPLE_STATES = ['solution', 'solid', 'bicelle', 'emulsion', 'fiber',
                 'filamentous virus', 'gel', 'liposome', 'membrane', 'micelle',
                 'lyophilized powder', 'oriented membrane film', 'fibrous protein',
                 'polycrystalline powder', 'reverse micelle', 'single crystal']

SHIFT_REF_COMPOUNDS = ['DSS','TSP','TMS','HDO','CHCl3','DMSO','DMSO-d5',
                       'dioxane','CH3OH','acetone','acetone-d5','TFE','TFA',
                       'H3PO4','TMP','PCr','liquid NH3']                 

STANDARD_UNITS = ['ppm','ppt','ppb']

NMR_PROBE_TYPES = ['liquid', 'solid', 'nano', 'flow', 'MAS']

GENERATION_TYPES = ['denovo', 'refinement']

SAMPLE_UNITS = ['kg', 'L', 'number']

SOLVENT_SYSTEMS = ['90% H2O/10% D2O',
                   '93% H2O/7% D2O',
                   '95% H2O/5% D2O',
                   '50% H2O/50% D2O',
                   '100% D2O',
                   'acetone',
                   'chloroform',
                   'DMSO',
                   'ethanol/water',
                   'methanol',
                   'trifluoroacetic acid',
                   'trifluoroethanol/water']

SAMPLE_COMP_UNITS = ['kg/m3', 'M', 'm3/m3', 'mol/mol', 'kg/kg']

# ionic strength not in API, but is useful to have it here

SAMPLE_CONDITION_TYPES = ['temperature', 'pressure', 'pH', 'ionic strength', 'spin-lock field',
                          'spin-lock offset', 'initial temperature', 'final temperature',
                          'delay time', 'mixing time']

POLYMER_PHYSICAL_STATES = ['native','denatured','molten globule','unfolded']

DB_NAMES = ['UniProt','GenBank','NDB','EMBL','PDB','BMRB','PIR','SWS']

NICE_GREEN = '#B0FFB0'

CITATION_STATUS = ['published','in press','in preparation']

JOURNALS = ['Arch. Biochem. Biophys.', 'Biochemistry',
            'Biochem. Biophys. Res. Commun.', 'Biochem. J.',
            'Biochim. Biophys. Acta', 'Biomol. NMR Assignments',
            'Biophys. J.', 'Biopolymers', 'Cell', 'EMBO J.' ,
            'Eur. J. Biochem.', 'FEBS Lett.', 'Inorg. Chem.',
            'Int. J. Pept. Protein Res.', 'J. Am. Chem. Soc.',
            'J. Biochem.', 'J. Biol. Chem.', 'J. Biol. Inorg. Chem.',
            'J. Biomol. NMR', 'J. Biomol. Struct. Dyn.',
            'J. Inorg. Biochem.', 'J. Magn. Reson.', 'J. Mol. Biol.',
            'J. Pept. Res.', 'J. Protein Chem.', 'Mol. Cell', 'Nature',
            'Nat. Struct. Biol.', 'Nucleic Acids Res.',
            'Proc. Natl. Acad. Sci. U.S.A.', 'Protein Eng.',
            'Proteins', 'Protein Sci.', 'Proteins: Struct. Funct. Genet.',
            'RNA', 'Science', 'Structure', 'Structure (Cambridge, MA, U.S.)']

CITATION_DEPOSIT = ['primary','other','no']

PRODUCTION_METHODS = ['recombinant technology', 'purified from the natural source',
                      'chemical synthesis', 'cell free synthesis', 'enzymatic semisynthesis',
                      'reverse transcriptase', 'obtained from a vendor']

VECTOR_TYPES = ['plasmid', 'baculovirus', 'other virus', 'cosmid']

ORGANISM_TYPES = ['organism', 'virus', 'plasmid', 'cosmid', 'no natural source']

"""
TOP 100:

HUMAN E   9606: N=Homo sapiens                C=Human
ECOLX B    562: N=Escherichia coli
MOUSE E  10090: N=Mus musculus                C=Mouse
YEAST E   4932: N=Saccharomyces cerevisiae                C=Baker's yeast
BOVIN E   9913: N=Bos taurus                C=Bovine
THET8 B 300852: N=Thermus thermophilus (strain HB8 / ATCC 27634 / DSM 579)
HALMA A   2238: N=Haloarcula marismortui
RAT   E  10116: N=Rattus norvegicus                C=Rat
THETH B    274: N=Thermus thermophilus
ECOLI B  83333: N=Escherichia coli (strain K12)

CHICK E   9031: N=Gallus gallus                C=Chicken
BACSU B   1423: N=Bacillus subtilis
BPT4  V  10665: N=Enterobacteria phage T4
Human immunodeficiency virus 1 (V 11676)
PIG   E   9823: N=Sus scrofa                C=Pig
THEMA B   2336: N=Thermotoga maritima
STAAU B   1280: N=Staphylococcus aureus
MYCTU B   1773: N=Mycobacterium tuberculosis
ARATH E   3702: N=Arabidopsis thaliana                C=Mouse-ear cress
PSEAE B    287: N=Pseudomonas aeruginosa

DEIRA B   1299: N=Deinococcus radiodurans
DROME E   7227: N=Drosophila melanogaster                C=Fruit fly
Geobacillus stearothermophilus (B 1422 Bacillus stearothermophilus)
RABIT E   9986: N=Oryctolagus cuniculus                C=Rabbit
Salmonella enterica subsp. enterica serovar Typhimurium (B 90371)
PSEPU B    303: N=Pseudomonas putida
RHOSH B   1063: N=Rhodobacter sphaeroides
Mycobacterium tuberculosis H37Rv (B 83332)
Thermus thermophilus HB27 (B 262724)
Pyrococcus horikoshii OT3 (A 70601)

METJA A   2190: N=Methanocaldococcus jannaschii
AQUAE B  63363: N=Aquifex aeolicus
ARCFU A   2234: N=Archaeoglobus fulgidus
SULSO A   2287: N=Sulfolobus solfataricus
PYRFU A   2261: N=Pyrococcus furiosus
PHYCA E   9755: N=Physeter catodon                C=Sperm whale
XENLA E   8355: N=Xenopus laevis                C=African clawed frog
HAEIN B    727: N=Haemophilus influenzae
HORSE E   9796: N=Equus caballus                C=Horse
Canis lupus familiaris (E 9615 Dog)

PYRHO A  53953: N=Pyrococcus horikoshii
PLAFA E   5833: N=Plasmodium falciparum
AEQVI E   6100: N=Aequorea victoria                C=Jellyfish
SPIOL E   3562: N=Spinacia oleracea                C=Spinach
Pseudomonas aeruginosa PAO1 (B 208964)
BACAM B   1390: N=Bacillus amyloliquefaciens
SCHPO E   4896: N=Schizosaccharomyces pombe                C=Fission yeast
HELPY B    210: N=Helicobacter pylori
VIBCH B    666: N=Vibrio cholerae
STRPN B   1313: N=Streptococcus pneumoniae

HRV14 V  12131: N=Human rhinovirus 14                C=HRV-14
Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 (B 99287)
Thermotoga maritima MSB8 (B 243274)
ALCFA B    511: N=Alcaligenes faecalis
Sulfolobus solfataricus P2 (A 273057)
STRAV B   1895: N=Streptomyces avidinii
CAEEL E   6239: N=Caenorhabditis elegans
BACCE B   1396: N=Bacillus cereus
AZOVI B    354: N=Azotobacter vinelandii
TRYBB E   5702: N=Trypanosoma brucei brucei

THEAQ B    271: N=Thermus aquaticus
Bacteroides thetaiotaomicron VPI-5482 (B 226186)
PARDE B    266: N=Paracoccus denitrificans
SOYBN E   3847: N=Glycine max                C=Soybean
THEAC A   2303: N=Thermoplasma acidophilum
Agrobacterium tumefaciens str. C58 (B 176299)
TRYCR E   5693: N=Trypanosoma cruzi
Methanothermobacter thermautotrophicus (A 145262)
SYNEL B  32046: N=Synechococcus elongatus
Influenza A virus (V 11320)

CVHSA V 227859: N=Human SARS coronavirus                C=SARS-CoV
Plasmodium falciparum 3D7 (E 36329)
Synechocystis sp. PCC 6803 (B 1148)
CLOBO B   1491: N=Clostridium botulinum
MAIZE E   4577: N=Zea mays                C=Maize
BACAN B   1392: N=Bacillus anthracis
KLEAE B  28451: N=Klebsiella aerogenes
DICDI E  44689: N=Dictyostelium discoideum                C=Slime mold
Halobacterium salinarum (A 2242)
HIRME E   6421: N=Hirudo medicinalis                C=Medicinal leech

PEA   E   3888: N=Pisum sativum                C=Garden pea
CHLRE E   3055: N=Chlamydomonas reinhardtii
PSEFL B    294: N=Pseudomonas fluorescens
HIV-1 M:B_HXB2R (V 11706 Human immunodeficiency virus type 1 (HXB2 ISOLATE))
HORVU E   4513: N=Hordeum vulgare                C=Barley
ENTFA B   1351: N=Enterococcus faecalis
TORCA E   7787: N=Torpedo californica                C=Pacific electric ray
RATRT E  10117: N=Rattus rattus                C=Black rat
SHEEP E   9940: N=Ovis aries                C=Sheep
Lactococcus lactis (B 1358)

Aquifex aeolicus VF5 (B 224324)
Human poliovirus 1 Mahoney (V 12081)
BACTH B   1427: N=Bacillus thermoproteolyticus
KLEPN B    573: N=Klebsiella pneumoniae
CLOTM B   1515: N=Clostridium thermocellum
ASPOR E   5062: N=Aspergillus oryzae
LACCA B   1582: N=Lactobacillus casei
METCA B    414: N=Methylococcus capsulatus
STRCL B   1901: N=Streptomyces clavuligerus
SERMA B    615: N=Serratia marcescens

NEXT ~150:

Archaeoglobus fulgidus DSM 4304  (A 224325)
THELA E   5541: N=Thermomyces lanuginosus
CANEN E   3823: N=Canavalia ensiformis                C=Jack bean
Streptomyces coelicolor A3(2) (B 100226)
STRCO B   1902: N=Streptomyces coelicolor
EMENI E 162425: N=Emericella nidulans
TRYBB E   5702: N=Trypanosoma brucei brucei
ZYMMO B    542: N=Zymomonas mobilis
ECO57 B  83334: N=Escherichia coli O157:H7
MELGA E   9103: N=Meleagris gallopavo                C=Common turkey
Hepatitis C virus (V 11103) 
Streptococcus pneumoniae TIGR4 (B 170187)
DANRE E   7955: N=Danio rerio                C=Zebrafish
Nostoc sp. PCC 7120 (B 103690)
Hypocrea jecorina (E 51453)
Blastochloris viridis (B 1079)
CLOPE B   1502: N=Clostridium perfringens
Streptomyces griseus (B 1911)
STRGR B   1911: N=Streptomyces griseus
Enterococcus faecalis V583 (B 226185)
SULTO A 111955: N=Sulfolobus tokodaii
RHOPR E  13249: N=Rhodnius prolixus                C=Triatomid bug
PYRAB A  29292: N=Pyrococcus abyssi
STRPY B   1314: N=Streptococcus pyogenes
synthetic construct (O 32630)
WHEAT E   4565: N=Triticum aestivum                C=Wheat
ARTGO B   1665: N=Arthrobacter globiformis
BACME B   1404: N=Bacillus megaterium
BACCI B   1397: N=Bacillus circulans
Helicobacter pylori 26695 (B 85962)
LAMBD V  10710: N=Enterobacteria phage lambda
STRLI B   1916: N=Streptomyces lividans
RHOCA B   1061: N=Rhodobacter capsulatus
Desulfovibrio vulgaris str. 'Miyazaki F' (B 883)
ERWCH B    556: N=Erwinia chrysanthemi
Bacillus halodurans C-125 (B  272558)
HALHA B   1053: N=Halorhodospira halophila
THEVL B  32053: N=Thermosynechococcus vulcanus
MASLA B  83541: N=Mastigocladus laminosus
Argopecten irradians (E 31199 Bay scallop)
RHIRD B    358: N=Rhizobium radiobacter                C=Agrobacterium tumefaciens
XANCP B    340: N=Xanthomonas campestris pv. campestris
CANAL E   5476: N=Candida albicans                C=Yeast
LEIMA E   5664: N=Leishmania major
BACHD B  86665: N=Bacillus halodurans
SHIFL B    623: N=Shigella flexneri
TOBAC E   4097: N=Nicotiana tabacum                C=Common tobacco
BACLI B   1402: N=Bacillus licheniformis
YERPE B    632: N=Yersinia pestis
BPT7  V  10760: N=Enterobacteria phage T7
Desulfovibrio vulgaris str. Hildenborough (B 882)
Methanothermobacter thermautotrophicus str. Delta H
THEEB B 197221: N=Thermosynechococcus elongatus (strain BP-1)
CORDI B   1717: N=Corynebacterium diphtheriae
FUSOX E   5507: N=Fusarium oxysporum                C=Panama disease fungus
CAMJE B    197: N=Campylobacter jejuni
Deinococcus radiodurans R1 (B 243230)
Engyodontium album (E 37998)
Acinetobacter sp. ADP1 (B 62977)
Achromobacter xylosoxidans (B 85698)
TOXGO E   5811: N=Toxoplasma gondii
Vaccinia virus (V  10245)
BURPS B  28450: N=Burkholderia pseudomallei
CORGL B   1718: N=Corynebacterium glutamicum
CAMDR E   9838: N=Camelus dromedarius                C=Dromedary
BPP22 V  10754: N=Enterobacteria phage P22
RHORU B   1085: N=Rhodospirillum rubrum
Daboia russellii pulchella (E 97228 Snake)
DESVU B    881: N=Desulfovibrio vulgaris
Nectria haematococca mpVI (E 70791)
unidentified (O 32644)
Salmonella enterica (B 28901)
RHOER B   1833: N=Rhodococcus erythropolis
PSEST B    316: N=Pseudomonas stutzeri
ENTCL B    550: N=Enterobacter cloacae
LYSEN B     69: N=Lysobacter enzymogenes
NEIME B    487: N=Neisseria meningitidis
Bacteroides fragilis NCTC 9343 (B 272559)
Escherichia coli BL21(DE3) (B 469008)
SULSH A   2286: N=Sulfolobus shibatae
YEREN B    630: N=Yersinia enterocolitica
Bacillus cereus ATCC 14579 (B 226900)
SCAIN E   6561: N=Scapharca inaequivalvis                C=Ark clam
Nostoc sp. PCC 7119 (B 1168)
APLCA E   6500: N=Aplysia californica                C=California sea hare
AERPE A  56636: N=Aeropyrum pernix
ASPNG E   5061: N=Aspergillus niger
BACTN B    818: N=Bacteroides thetaiotaomicron
Burkholderia pseudomallei 1710b (B 320372)
STRAU B   1894: N=Streptomyces aureofaciens
RHOPA B   1076: N=Rhodopseudomonas palustris
RICCO E   3988: N=Ricinus communis                C=Castor bean
Shewanella oneidensis MR-1 (B 211586)
STRMU B   1309: N=Streptococcus mutans
MYCSM B   1772: N=Mycobacterium smegmatis
Streptococcus pneumoniae R6 (B 171101)
CRYPA E   5116: N=Cryphonectria parasitica                C=Chesnut blight fungus
PSEME B    300: N=Pseudomonas mendocina
DESDE B    876: N=Desulfovibrio desulfuricans
Escherichia coli str. K-12 substr. W3110 (B 316407)
Enterobacterio phage MS2 (V 12022)
Cellvibrio japonicus (B 155077)
CLOPA B   1501: N=Clostridium pasteurianum
SULAC A   2285: N=Sulfolobus acidocaldarius
BRAJA B    375: N=Bradyrhizobium japonicum
PSESP B    306: N=Pseudomonas sp.
Mycobacterium smegmatis str. MC2 155 (B 246196)
Thermococcus kodakarensis KOD1 (A 69014)
Synechococcus elongatus PCC 7942 (B 1140) 
Rhodococcus jostii RHA1 (B 101510)
STEMA B  40324: N=Stenotrophomonas maltophilia                C=Pseudomonas maltophilia
VISAL E   3972: N=Viscum album                C=European mistletoe
BORPE B    520: N=Bordetella pertussis
BPPH2 V  10756: N=Bacillus phage phi29
LACPL B   1590: N=Lactobacillus plantarum
Bacillus anthracis str. Ames (B 198094)
Sulfolobus tokodaii str. 7 (A 273063)
Sinorhizobium meliloti (B 382) 
BURCE B    292: N=Burkholderia cepacia
THEVU B   2026: N=Thermoactinomyces vulgaris
MEDSA E   3879: N=Medicago sativa                C=Alfalfa
NEUCR E   5141: N=Neurospora crassa
Listeria monocytogenes EGD-e (B 169963)
PROFR B   1752: N=Propionibacterium freudenreichii subsp. shermanii
Pyrococcus furiosus DSM 3638 (A 186497)
MESMA E  34649: N=Mesobuthus martensii                C=Manchurian scorpion
PHACH E   5306: N=Phanerochaete chrysosporium                C=White-rot fungus
Escherichia coli B (B 37762)
NEIGO B    485: N=Neisseria gonorrhoeae
BORBU B    139: N=Borrelia burgdorferi                C=Lyme disease spirochete
Staphylococcus aureus subsp. aureus Mu50 (B 158878)
ARMRU E   3704: N=Armoracia rusticana                C=Horseradish
Burkholderia xenovorans LB400 (B 266265)
PYRAE A  13773: N=Pyrobaculum aerophilum
Thermosynechococcus elongatus BP-1 (B 197221)
LISMO B   1639: N=Listeria monocytogenes
PLAVI E   5855: N=Plasmodium vivax
ORYSA E   4530: N=Oryza sativa                C=Rice
HEVBR E   3981: N=Hevea brasiliensis                C=Para rubber tree
LEIME E   5665: N=Leishmania mexicana
NITEU B    915: N=Nitrosomonas europaea
CAPHI E   9925: N=Capra hircus                C=Goat
ARTIN E   3490: N=Artocarpus integer                C=Jack fruit
Porphyromonas gingivalis W83 (B 242619)
Aeropyrum pernix K1 (A 272557)
MACAM E   8199: N=Macrozoarces americanus                C=Ocean pout
MLVMO V  11801: N=Moloney murine leukemia virus                C=MoMLV
DESGI B    879: N=Desulfovibrio gigas
Methanocaldococcus jannaschii DSM 2661 (A 243232)
Poliovirus type 3 (strains P3/LEON/37 AND P3/LEON 12A[1]B) (V 12088 Human Poliovirus type 3)
Sterkiella nova (E 200597)

"""

ORGANISMS = {'Eukarya A-E':  [('Aequorea victoria', 'Jellyfish'),
                              ('Anopheles gambiae','Mosquito'),
                              ('Aplysia californica', 'California sea hare'),
                              ('Arabidopsis thaliana', 'Mouse-ear cress'),
                              ('Armoracia rusticana', 'Horseradish'),
                              ('Artocarpus integer', 'Jack fruit'),
                              ('Aspergillus niger', 'Aspergillus niger'),
                              ('Aspergillus oryzae', 'Aspergillus oryzae'),
                              ('Bos taurus', 'Cow'),
                              ('Caenorhabditis elegans', 'Caenorhabditis elegans'),
                              ('Camelus dromedarius', 'Dromedary'),
                              ('Canavalia ensiformis', 'Jack bean'),
                              ('Candida albicans', 'Yeast'),
                              ('Canis lupus familiaris', 'Dog'),
                              ('Capra hircus', 'Goat'),
                              ('Chlamydomonas reinhardtii', 'Chlamydomonas reinhardtii'),
                              ('Cryphonectria parasitica', 'Chesnut blight fungus'),
                              ('Daboia russellii pulchella', 'Snake'),
                              ('Danio rerio', 'Zebrafish'),
                              ('Dictyostelium discoideum', 'Slime mold'),
                              ('Drosophila melanogaster', 'Fruit fly'),
                              ('Emericella nidulans', 'Emericella nidulans'),
                              ('Engyodontium album', 'Engyodontium album'),
                              ('Equus caballus', 'Horse'),
                              ],
             'Eukarya F-J':  [('Fusarium oxysporum', 'Panama disease fungus'),
                              ('Gallus gallus', 'Chicken'),
                              ('Glycine max', 'Soybean'),
                              ('Hevea brasiliensis', 'Para rubber tree'),
                              ('Hirudo medicinalis', 'Medicinal leech'),
                              ('Homo sapiens', 'Human'),
                              ('Hordeum vulgare', 'Barley'),
                              ('Hypocrea jecorina', 'Hypocrea jecorina'),
                              ],
             'Eukarya K-O':  [('Leishmania major', 'Leishmania major'),
                              ('Leishmania mexicana', 'Leishmania mexicana'),
                              ('Macrozoarces americanus', 'Ocean pout'),
                              ('Medicago sativa', 'Alfalfa'),
                              ('Meleagris gallopavo', 'Common turkey'),
                              ('Mesobuthus martensii', 'Manchurian scorpion'),
                              ('Mus musculus', 'Mouse'),
                              ('Nectria haematococca mpVI', 'Nectria haematococca mpVI'),
                              ('Neurospora crassa', 'Neurospora crassa'),
                              ('Nicotiana tabacum', 'Common tobacco'),
                              ('Oryctolagus cuniculus', 'Rabbit'),
                              ('Oryza sativa', 'Rice'),
                              ('Ovis aries', 'Sheep'),
                              ],
             'Eukarya P-T':  [('Phanerochaete chrysosporium', 'White-rot fungus'),
                              ('Physeter catodon', 'Sperm whale'),
                              ('Pisum sativum', 'Garden pea'),
                              ('Plasmodium falciparum', 'Plasmodium falciparum'),
                              ('Plasmodium vivax', 'Plasmodium vivax'),
                              ('Plasmodium falciparum 3D7', 'Plasmodium falciparum 3D7'),
                              ('Rattus norvegicus', 'Rat'),
                              ('Rattus rattus', 'Black rat'),
                              ('Rhodnius prolixus', 'Triatomid bug'),
                              ('Ricinus communis', 'Castor bean'),
                              ('Saccharomyces cerevisiae', "Baker's yeast"),
                              ('Scapharca inaequivalvis', 'Ark clam'),
                              ('Schizosaccharomyces pombe', 'Fission yeast'),
                              ('Spinacia oleracea', 'Spinach'),
                              ('Sus scrofa', 'Pig'),
                              ('Sterkiella nova', 'Sterkiella nova'),
                              ('Thermomyces lanuginosus', 'Thermomyces lanuginosus'),
                              ('Torpedo californica', 'Pacific electric ray'),
                              ('Toxoplasma gondii', 'Toxoplasma gondii'),
                              ('Trypanosoma brucei brucei', 'Trypanosoma brucei brucei'),
                              ('Trypanosoma cruzi', 'Trypanosoma cruzi'),
                              ('Triticum aestivum', 'Wheat'),
                              ('Trypanosoma brucei brucei', 'Trypanosoma brucei brucei'),
                              ],
             'Eukarya V-Z':  [('Viscum album', 'European mistletoe'),
                              ('Xenopus laevis', 'African clawed frog'),
                              ('Zea mays', 'Maize'),
                              ],
             'Bacteria A-B': [('Achromobacter xylosoxidans', 'Achromobacter xylosoxidans'),
                              ('Acinetobacter sp. ADP1', 'Acinetobacter sp. ADP1'),
                              ('Agrobacterium tumefaciens', 'Agrobacterium tumefaciens'),
                              ('Agrobacterium tumefaciens str. C58', 'Agrobacterium tumefaciens str. C58'),
                              ('Alcaligenes faecalis', 'Alcaligenes faecalis'),
                              ('Aquifex aeolicus', 'Aquifex aeolicus'),
                              ('Aquifex aeolicus VF5', 'Aquifex aeolicus VF5'),
                              ('Arthrobacter globiformis', 'Arthrobacter globiformis'),
                              ('Azotobacter vinelandii', 'Azotobacter vinelandii'),
                              ('Bacillus amyloliquefaciens', 'Bacillus amyloliquefaciens'),
                              ('Bacillus anthracis', 'Bacillus anthracis'),
                              ('Bacillus anthracis str. Ames', 'Bacillus anthracis str. Ames'),
                              ('Bacillus cereus', 'Bacillus cereus'),
                              ('Bacillus cereus ATCC 14579', 'Bacillus cereus ATCC 14579'),
                              ('Bacillus subtilis', 'Bacillus subtilis'),
                              ('Bacillus circulans', 'Bacillus circulans'),
                              ('Bacillus halodurans', 'Bacillus halodurans'),
                              ('Bacillus licheniformis', 'Bacillus licheniformis'),
                              ('Bacillus megaterium', 'Bacillus megaterium'),
                              ('Bacillus thermoproteolyticus', 'Bacillus thermoproteolyticus'),
                              ('Bacteroides thetaiotaomicron VPI-5482', 'Bacteroides thetaiotaomicron VPI-5482'),
                              ('Bacteroides fragilis NCTC 9343', 'Bacteroides fragilis NCTC 9343'),
                              ('Bacteroides thetaiotaomicron', 'Bacteroides thetaiotaomicron'),
                              ('Blastochloris viridis', 'Blastochloris viridis'),
                              ('Bordetella pertussis', 'Bordetella pertussis'),
                              ('Borrelia burgdorferi', 'Lyme disease spirochete'),
                              ('Bradyrhizobium japonicum', 'Bradyrhizobium japonicum'),
                              ('Burkholderia cepacia', 'Burkholderia cepacia'),
                              ('Burkholderia pseudomallei', 'Burkholderia pseudomallei'),
                              ('Burkholderia pseudomallei 1710b', 'Burkholderia pseudomallei 1710b'),
                              ('Burkholderia xenovorans LB400', 'Burkholderia xenovorans LB400'),
                              ],
             'Bacteria C-E': [('Campylobacter jejuni', 'Campylobacter jejuni'),
                              ('Cellvibrio japonicus', 'Cellvibrio japonicus'),
                              ('Clostridium botulinum', 'Clostridium botulinum'),
                              ('Clostridium pasteurianum', 'Clostridium pasteurianum'),
                              ('Clostridium perfringens', 'Clostridium perfringens'),
                              ('Clostridium thermocellum', 'Clostridium thermocellum'),
                              ('Corynebacterium diphtheriae', 'Corynebacterium diphtheriae'),
                              ('Corynebacterium glutamicum', 'Corynebacterium glutamicum'),
                              ('Deinococcus radiodurans', 'Deinococcus radiodurans'),
                              ('Deinococcus radiodurans R1', 'Deinococcus radiodurans R1'),
                              ('Desulfovibrio desulfuricans', 'Desulfovibrio desulfuricans'),
                              ('Desulfovibrio gigas', 'Desulfovibrio gigas'),
                              ('Desulfovibrio vulgaris', 'Desulfovibrio vulgaris'),
                              ("Desulfovibrio vulgaris str. 'Miyazaki F'", "Desulfovibrio vulgaris str. 'Miyazaki F'"),
                              ('Desulfovibrio vulgaris str. Hildenborough', 'Desulfovibrio vulgaris str. Hildenborough'),
                              ('Enterococcus faecalis', 'Enterococcus faecalis'),
                              ('Enterobacter cloacae', 'Enterobacter cloacae'),
                              ('Enterococcus faecalis V583', 'Enterococcus faecalis V583'),
                              ('Erwinia chrysanthemi', 'Erwinia chrysanthemi'),
                              ('Erwinia tasmaniensis', 'Erwinia tasmaniensis'),
                              ('Erythrobacter litoralis', 'Erythrobacter litoralis'),
                              ('Escherichia coli', 'Escherichia coli'),
                              ('Escherichia coli B', 'Escherichia coli B'),
                              ('Escherichia coli BL21(DE3)', 'Escherichia coli BL21(DE3)'),
                              ('Escherichia coli O157:H7', 'Escherichia coli O157:H7'),
                              ('Escherichia coli (strain K12)', 'Escherichia coli (strain K12)'),
                              ('Escherichia coli str. K-12 substr. W3110', 'Escherichia coli str. K-12 substr. W3110'),
                              ('Escherichia fergusonii', 'Escherichia fergusonii'),
                              ('Exiguobacterium sibiricum', 'Exiguobacterium sibiricum'),
                              ],
             'Bacteria F-L': [('Fervidobacterium nodosum', 'Fervidobacterium nodosum'),
                              ('Finegoldia magna', 'Finegoldia magna'),
                              ('Flavobacterium johnsoniae', 'Flavobacterium johnsoniae'),
                              ('Flavobacterium psychrophilu', 'Flavobacterium psychrophilu'),
                              ('Francisella novicida', 'Francisella novicida'),
                              ('Francisella philomiragia', 'Francisella philomiragia'),
                              ('Francisella tularensis', 'Francisella tularensis'),
                              ('Frankia alni', 'Frankia alni'),
                              ('Fusobacterium nucleatum', 'Fusobacterium nucleatum'),
                              ('Geobacillus kaustophilus', 'Geobacillus kaustophilus'),
                              ('Geobacillus thermodenitrificans', 'Geobacillus thermodenitrificans'),
                              ('Geobacter bemidjiensis', 'Geobacter bemidjiensis'),
                              ('Haemophilus influenzae', 'Haemophilus influenzae'),
                              ('Halorhodospira halophila', 'Halorhodospira halophila'),
                              ('Helicobacter pylori', 'Helicobacter pylori'),
                              ('Helicobacter pylori 26695', 'Helicobacter pylori 26695'),
                              ('Klebsiella aerogenes', 'Klebsiella aerogenes'),
                              ('Klebsiella pneumoniae', 'Klebsiella pneumoniae'),
                              ('Lactobacillus casei', 'Lactobacillus casei'),
                              ('Lactococcus lactis', 'Lactococcus lactis'),
                              ('Lactobacillus plantarum', 'Lactobacillus plantarum'),
                              ('Listeria monocytogenes', 'Listeria monocytogenes'),
                              ('Listeria monocytogenes EGD-e', 'Listeria monocytogenes EGD-e'),
                              ('Lysobacter enzymogenes', 'Lysobacter enzymogenes'),
                              ],
             'Bacteria M-P': [('Mastigocladus laminosus', 'Mastigocladus laminosus'),
                              ('Mycobacterium smegmatis', 'Mycobacterium smegmatis'),
                              ('Mycobacterium smegmatis str. MC2 155', 'Mycobacterium smegmatis str. MC2 155'),
                              ('Methylococcus capsulatus', 'Methylococcus capsulatus'),
                              ('Mycobacterium tuberculosis', 'Mycobacterium tuberculosis'),
                              ('Mycobacterium tuberculosis H37Rv', 'Mycobacterium tuberculosis H37Rv'),
                              ('Neisseria gonorrhoeae', 'Neisseria gonorrhoeae'),
                              ('Neisseria meningitidis', 'Neisseria meningitidis'),
                              ('Nitrosomonas europaea', 'Nitrosomonas europaea'),
                              ('Nostoc sp. PCC 7119', 'Nostoc sp. PCC 7119'),
                              ('Nostoc sp. PCC 7120', 'Nostoc sp. PCC 7120'),
                              ('Paracoccus denitrificans', 'Paracoccus denitrificans'),
                              ('Porphyromonas gingivalis W83', 'Porphyromonas gingivalis W83'),
                              ('Propionibacterium freudenreichii subsp. shermanii', 'Propionibacterium freudenreichii subsp. shermanii'),
                              ('Pseudomonas aeruginosa', 'Pseudomonas aeruginosa'),
                              ('Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa PAO1'),
                              ('Pseudomonas fluorescens', 'Pseudomonas fluorescens'),
                              ('Pseudomonas putida', 'Pseudomonas putida'),
                              ('Pseudomonas mendocina', 'Pseudomonas mendocina'),
                              ('Pseudomonas sp.', 'Pseudomonas sp.'),
                              ('Pseudomonas stutzeri', 'Pseudomonas stutzeri'),
                              ],
             'Bacteria Q-S': [('Rhizobium radiobacter', 'Agrobacterium tumefaciens'),
                              ('Rhodobacter capsulatus', 'Rhodobacter capsulatus'),
                              ('Rhodobacter sphaeroides', 'Rhodobacter sphaeroides'),
                              ('Rhodococcus erythropolis', 'Rhodococcus erythropolis'),
                              ('Rhodococcus jostii RHA1', 'Rhodococcus jostii RHA1'),
                              ('Rhodopseudomonas palustris', 'Rhodopseudomonas palustris'),
                              ('Rhodospirillum rubrum', 'Rhodospirillum rubrum'),
                              ('Salmonella enterica', 'Salmonella enterica'),
                              ('Salmonella enterica subsp. enterica serovar Typhimurium', 'Salmonella enterica subsp. enterica serovar Typhimurium'),
                              ('Salmonella enterica subsp. enterica serovar Typhimurium str. LT2', 'Salmonella enterica subsp. enterica serovar Typhimurium str. LT2'),
                              ('Serratia marcescens', 'Serratia marcescens'),
                              ('Shewanella oneidensis MR-1', 'Shewanella oneidensis MR-1'),
                              ('Shigella flexneri', 'Shigella flexneri'),
                              ('Staphylococcus aureus', 'Staphylococcus aureus'),
                              ('Staphylococcus aureus subsp. aureus Mu50', 'Staphylococcus aureus subsp. aureus Mu50'),
                              ('Stenotrophomonas maltophilia', 'Pseudomonas maltophilia'),
                              ('Streptococcus mutans', 'Streptococcus mutans'),
                              ('Streptococcus pneumoniae', 'Streptococcus pneumoniae'),
                              ('Streptococcus pneumoniae R6', 'Streptococcus pneumoniae R6'),
                              ('Streptococcus pneumoniae TIGR4', 'Streptococcus pneumoniae TIGR4'),
                              ('Streptococcus pyogenes', 'Streptococcus pyogenes'),
                              ('Streptomyces aureofaciens', 'Streptomyces aureofaciens'),
                              ('Streptomyces avidinii', 'Streptomyces avidinii'),
                              ('Streptomyces clavuligerus', 'Streptomyces clavuligerus'),
                              ('Streptomyces coelicolor', 'Streptomyces coelicolor'),
                              ('Streptomyces coelicolor A3(2)', 'Streptomyces coelicolor A3(2)'),
                              ('Streptomyces griseus', 'Streptomyces griseus'),
                              ('Streptomyces lividans', 'Streptomyces lividans'),
                              ('Synechococcus elongatus', 'Synechococcus elongatus'),
                              ('Synechocystis sp. PCC 6803', 'Synechocystis sp. PCC 6803'),
                              ],
             'Bacteria T-Z': [('Thermoactinomyces vulgaris', 'Thermoactinomyces vulgaris'),
                              ('Thermosynechococcus elongatus (strain BP-1)', 'Thermosynechococcus elongatus (strain BP-1)'),
                              ('Thermosynechococcus elongatus BP-1', 'Thermosynechococcus elongatus BP-1'),
                              ('Thermosynechococcus vulcanus', 'Thermosynechococcus vulcanus'),
                              ('Thermotoga maritima', 'Thermotoga maritima'),
                              ('Thermotoga maritima MSB8', 'Thermotoga maritima MSB8'),
                              ('Thermus aquaticus', 'Thermus aquaticus'),
                              ('Thermus thermophilus', 'Thermus thermophilus'),
                              ('Thermus thermophilus HB27', 'Thermus thermophilus HB27'),
                              ('Thermus thermophilus (strain HB8 / ATCC 27634 / DSM 579)', 'Thermus thermophilus (strain HB8 / ATCC 27634 / DSM 579)'),
                              ('Vibrio cholerae', 'Vibrio cholerae'),
                              ('Xanthomonas campestris pv. campestris', 'Xanthomonas campestris pv. campestris'),
                              ('Yersinia enterocolitica', 'Yersinia enterocolitica'),
                              ('Yersinia pestis', 'Yersinia pestis'),
                              ('Zymomonas mobilis', 'Zymomonas mobilis'),
                          ],
             'Archaea A-M':  [('Aeropyrum pernix', 'Aeropyrum pernix'),
                              ('Aeropyrum pernix K1', 'Aeropyrum pernix K1'),
                              ('Archaeoglobus fulgidus', 'Archaeoglobus fulgidus'),
                              ('Archaeoglobus fulgidus DSM 4304 ', 'Archaeoglobus fulgidus DSM 4304 '),
                              ('Caldivirga maquilingensis', 'Caldivirga maquilingensis'),
                              ('Candidatus Korarchaeum', 'Candidatus Korarchaeum'),
                              ('Candidatus Methanoregula', 'Candidatus Methanoregula'),
                              ('Candidatus Methanosphaerula', 'Candidatus Methanosphaerula'),
                              ('Desulfurococcus kamchatkensis', 'Desulfurococcus kamchatkensis'),
                              ('Haloarcula marismortui', 'Haloarcula marismortui'),
                              ('Halobacterium salinarum', 'Halobacterium salinarum'),
                              ('Halobacterium sp.', 'Halobacterium sp.'),
                              ('Haloquadratum walsbyi', 'Haloquadratum walsbyi'),
                              ('Halorubrum lacusprofundi', 'Halorubrum lacusprofundi'),
                              ('Hyperthermus butylicus', 'Hyperthermus butylicus'),
                              ('Ignicoccus hospitalis', 'Ignicoccus hospitalis'),
                              ('Metallosphaera sedula', 'Metallosphaera sedula'),
                              ('Methanobrevibacter smithii', 'Methanobrevibacter smithii'),
                              ('Methanocaldococcus jannaschii', 'Methanocaldococcus jannaschii'),
                              ('Methanocaldococcus jannaschii DSM 2661', 'Methanocaldococcus jannaschii DSM 2661'),
                              ('Methanococcoides burtonii', 'Methanococcoides burtonii'),
                              ('Methanococcus aeolicus', 'Methanococcus aeolicus'),
                              ('Methanococcus vannielii', 'Methanococcus vannielii'),
                              ('Methanocorpusculum labreanum', 'Methanocorpusculum labreanum'),
                              ('Methanoculleus marisnigri', 'Methanoculleus marisnigri'),
                              ('Methanopyrus kandleri', 'Methanopyrus kandleri'),
                              ('Methanosaeta thermophila', 'Methanosaeta thermophila'),
                              ('Methanosarcina acetivorans', 'Methanosarcina acetivorans'),
                              ('Methanosarcina barkeri', 'Methanosarcina barkeri'),
                              ('Methanosarcina mazei', 'Methanosarcina mazei'),
                              ('Methanosphaera stadtmanae', 'Methanosphaera stadtmanae'),
                              ('Methanospirillum hungatei', 'Methanospirillum hungatei'),
                              ('Methanothermobacter thermautotrophicus', 'Methanothermobacter thermautotrophicus'),
                              ],
             'Archaea N-Z':  [('Nanoarchaeum equitans', 'Nanoarchaeum equitans'),
                              ('Natronomonas pharaonis', 'Natronomonas pharaonis'),
                              ('Nitrosopumilus maritimus', 'Nitrosopumilus maritimus'),
                              ('Picrophilus torridus', 'Picrophilus torridus'),
                              ('Pyrobaculum arsenaticum', 'Pyrobaculum arsenaticum'),
                              ('Pyrobaculum calidifontis', 'Pyrobaculum calidifontis'),
                              ('Pyrobaculum islandicum', 'Pyrobaculum islandicum'),
                              ('Pyrobaculum aerophilum', 'Pyrobaculum aerophilum'),
                              ('Pyrococcus abyssi', 'Pyrococcus abyssi'),
                              ('Pyrococcus furiosus', 'Pyrococcus furiosus'),
                              ('Pyrococcus furiosus DSM 3638', 'Pyrococcus furiosus DSM 3638'),
                              ('Pyrococcus horikoshii', 'Pyrococcus horikoshii'),
                              ('Pyrococcus horikoshii OT3', 'Pyrococcus horikoshii OT3'),
                              ('Sulfolobus acidocaldarius', 'Sulfolobus acidocaldarius'),
                              ('Sulfolobus shibatae', 'Sulfolobus shibatae'),
                              ('Sulfolobus solfataricus', 'Sulfolobus solfataricus'),
                              ('Sulfolobus solfataricus P2', 'Sulfolobus solfataricus P2'),
                              ('Sulfolobus tokodaii', 'Sulfolobus tokodaii'),
                              ('Sulfolobus tokodaii str. 7', 'Sulfolobus tokodaii str. 7'),
                              ('Staphylothermus marinus', 'Staphylothermus marinus'),
                              ('Thermococcus kodakarensis', 'Thermococcus kodakarensis'),
                              ('Thermococcus kodakarensis KOD1', 'Thermococcus kodakarensis KOD1'),
                              ('Thermococcus onnurineus', 'Thermococcus onnurineus'),
                              ('Thermofilum pendens', 'Thermofilum pendens'),
                              ('Thermoplasma acidophilum', 'Thermoplasma acidophilum'),
                              ('Thermoplasma volcanium', 'Thermoplasma volcanium'),
                              ('Thermoproteus neutrophilus', 'Thermoproteus neutrophilus'),
                              ],
             'Other':        [('Synthetic construct', 'Synthetic construct'),
                              ('Unidentified', 'Unidentified'),
                              ],
             'Viruses':      [('Bacillus phage phi29', 'Bacillus phage phi29'),
                              ('Enterobacteria phage P22', 'Enterobacteria phage P22'),
                              ('Enterobacteria phage T4', 'Enterobacteria phage T4'),
                              ('Enterobacteria phage T7', 'Enterobacteria phage T7'),
                              ('Enterobacteria phage lambda', 'Enterobacteria phage lambda'),
                              ('Enterobacterio phage MS2', 'Enterobacterio phage MS2'),
                              ('Human immunodeficiency virus 1', 'Human immunodeficiency virus 1'),
                              ('Human poliovirus 1 Mahoney', 'Human poliovirus 1 Mahoney'),
                              ('Poliovirus type 3', 'Human poliovirus type 3'),
                              ('Human rhinovirus 14', 'HRV-14'),
                              ('Human SARS coronavirus', 'SARS-CoV'),                          
                              ('Influenza A virus', 'Influenza A virus'),
                              ('Moloney murine leukemia virus', 'MoMLV'),
                              ],
             }

EXP_ORGANISMS = ['Aspergillus niger', 'Cricetulus griseus', 'Desulfovibrio desulfuricans',
                 'Escherichia coli', 'Pichia pastoris', 'Saccharomyces cerevisiae',
                 'Spodoptera frugiperda', 'Synthetic', 'Unidentified',
                 'Wheat germ - cell free' 'E. coli - cell free']

SPEC_MANU_LIST = ['Varian', 'Bruker', 'JEOL', 'General Electric',
                  'Chemagnetics', 'Custom Made', 'Other']

SPEC_MODEL_LIST = [('Varian', 'UnityInova'),
                   ('Varian', 'MERCURYplus'),
                   ('Varian', 'InfinityPlus'),
                   #('Varian', 'Other'),
                   ('Bruker', 'AM'),
                   ('Bruker', 'AMX'),
                   ('Bruker', 'DMX'),
                   ('Bruker', 'Avance'),
                   #('Bruker', 'Other'),
                   ('JEOL', 'ECX'),
                   ('JEOL', 'ECA'),
                   #('JEOL', 'Other'),
                   #('General Electric', 'Other'),
                   #('Chemagnetics', 'Other'),
                   #('Custom Made', 'Other'),
                   #('Other', 'Other'),
                   ]

ISOTOPE_LABEL_LIST =  [('natural abundance', 'natural abundance'),
                       ('[U-15N]', '[U-100% 15N]'),
                       ('[U-15N]', '[U-99% 15N]'),
                       ('[U-15N]', '[U-98% 15N]'),
                       ('[U-15N]', '[U-95% 15N]'),
                       ('[U-15N]', '[U-90% 15N]'),
                       ('[U-15N]', '[U-15N]'),
                       ('[U-13C]', '[U-100% 13C]'),
                       ('[U-13C]', '[U-95% 13C]'),
                       ('[U-13C]', '[U-10% 13C]'),
                       ('[U-13C]', '[U-13C]'),
                       ('[U-13C; U-15N]', '[U-100% 13C; U-100% 15N]'),
                       ('[U-13C; U-15N]', '[U-99% 13C; U-99% 15N]'),
                       ('[U-13C; U-15N]', '[U-98% 13C; U-98% 15N]'),
                       ('[U-13C; U-15N]', '[U-95% 13C; U-95% 15N]'),
                       ('[U-13C; U-15N]', '[U-95% 13C; U-90% 15N]'),
                       ('[U-13C; U-15N]', '[U-10% 13C; U-100% 15N]'),
                       ('[U-13C; U-15N]', '[U-10% 13C; U-99% 15N]'),
                       ('[U-13C; U-15N]', '[U-13C; U-15N]'),
                       ('[U-13C; U-15N; U-2H]', '[U-100% 13C; U-100% 15N; U-80% 2H]'),
                       ('[U-13C; U-15N; U-2H]', '[U-95% 13C; U-95% 15N; U-95% 2H]'),
                       ('[U-13C; U-15N; U-2H]', '[U-13C; U-15N; U-2H]'),
                       ('[U-2H]', '[U-100% 2H]'),
                       ('[U-2H]', '[U-99% 2H]'),
                       ('[U-2H]', '[U-2H]'),
                       ('Residue-specific', '[U-13C; U-15N]-Ade'),
                       ('Residue-specific', '[U-13C; U-15N]-Cyt'),
                       ('Residue-specific', '[U-13C; U-15N]-Gua'),
                       ('Residue-specific', '[U-13C; U-15N]-Ura'),
                       ('Residue-specific', '[U-15N]-Leu'),
                       ('Atom-specific', '[95% 13CA]-Trp'),
                       ]

class EntryCompletionFrame(Frame):

  '''
ECI allows you to easily add an "Entry" object to your CCPN project. An "Entry" object contains all the information that you wish to deposit with your submission. You can also select chemical shift lists, peak lists, structural restraints, ensembles, etc. and all at the click of a button. In addition, you can add all the meta data that is required for submissions to the PDB and BMRB. This can be done securely on your desktop computer over the duration of your NMR project.
  '''

  def __init__(self, parent, basePopup, *args, **kw):

    self.waiting = False
    self.entryObjType = 'ShiftList'
    self.editObject = None

    self.project = basePopup.project

    if self.project:
      self.nmrProject = self.project.currentNmrProject
    else:
      self.nmrProject = None  

    self.entry = None

    self.submissionType = 'PDB'

    self.person = None
    self.author = None
    self.contactPerson = None
    self.group = None

    self.citation = None
    self.software = None

    self.shiftReference = None
    self.checkShiftList = None

    self.molSystem = None
    self.chain = None
    self.molecule = None
    self.align = None
    self.dbName = None
    self.accNum = None

    self.entryMolecule = None
    self.expSource = None
    self.natSource = None

    self.sample = None
    self.sampleComponent = None

    self.spectrometer = None
    self.probe = None
    self.experiment = None

    self.sampleConditionSet = None
    self.sampleCondition = None

    self.strucGen = None
    self.ensemble = None
    self.constSet = None
    self.constList = None

    self.registerNotify = basePopup.registerNotify
    self.unregisterNotify = basePopup.unregisterNotify

    self.classSubmitDict = {'PeakList':('PeakList','peakLists'),}

    for className in MEASUREMENT_LIST_CLASSES:
      self.classSubmitDict[className] = ('MeasurementList', 'measurementLists')

    for className in DERIVED_LIST_CLASSES:
      self.classSubmitDict[className] = ('DerivedDataList', 'derivedDataLists')

    Frame.__init__(self, parent, **kw)

    self.grid_columnconfigure(1, weight=1)
    self.grid_rowconfigure(1, weight=1)

    options = ['Main', 'People', 'References', 'NMR Data',
               'Shifts', 'Molecules', 'Sources', 'Samples',
               'NMR Conditions', 'Instruments', 'Experiments', 'Structures']

    tipTexts = ['Add a new "Entry" and check the completeness of it',
                'Add people, group and address information',
                'Add citations and software used',
                'Add chemical shift lists, peak lists and other NMR datasets',
                'Add chemical shift references and look for chemical shift outliers',
                'Review the molecules, chains and unusual chemical components in your molecular system',
                'Add experimental and biological sources',
                'Add sample, sample component and isotope labelling information',
                'Add experimental conditions used in your NMR experiments',
                'Add NMR spectrometers and probes',
                'Add NMR experiments and link data about samples, conditions, instruments, shift references and shift lists added in the previous tabs to the experiments',
                'Add a structure generation and link a structural ensemble and restraint sets to it'
                ]

    tabbedFrame = TabbedFrame(self, options=options, tipTexts=tipTexts, callback=self.selectTab)
    tabbedFrame.grid(row=1, column=0, columnspan=2, sticky='nsew')

    frameM, frameA, frameB, frameC, frameD, frameE, frameF, frameG, frameH, frameI, frameJ, frameK = tabbedFrame.frames

    self.tabbedFrame = tabbedFrame

    #
    # Main
    #

    i = 0

    frameM.grid_columnconfigure(1, weight=1)
    frameM.grid_rowconfigure(9, weight=1)

    label = Label(frameM, text='Deposition Entry')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Name of "Entry" being edited, add a new one with the button at the bottom left of this tab'

    self.entryPulldown = PulldownList(frameM, tipText=tipText, callback=self.changeEntry)
    self.entryPulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Molecular System')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Name of molecular system used in the study - review it in more detail in the "Molecules" tab'

    self.molSystemPulldown = PulldownList(frameM, tipText=tipText, callback=self.changeMolSystem)
    self.molSystemPulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Title')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Title that will be in the header of the PDB file'

    self.titleEntry = Entry(frameM, text='', width=32, tipText=tipText,
                            returnCallback=self.updateEntryTitle)
    self.titleEntry.grid(row=i, column=1, sticky='ew')
    self.titleEntry.bind('<Leave>', self.updateEntryTitle, '+') 
    i += 1

    label = Label(frameM, text='NMR Method Type')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Type of NMR data collected, choose from: solution, solid-state, theoretical or ?'

    nmrMethodTypes = ['solution','solid-state','theoretical','?']

    self.nmrMethodPulldown = PulldownList(frameM, texts=nmrMethodTypes, tipText=tipText,
                                          callback=self.setEntryType)
    self.nmrMethodPulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Submission Type')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Is the deposition for the PDB or BMRB?  The mandatory fields in the completeness report below will reflect the type of deposition to be done'

    submissionTypes = ['PDB', 'BMRB']

    self.subTypePulldown = PulldownList(frameM, texts=submissionTypes, tipText=tipText,
                                        callback=self.setSubmissionType)
    self.subTypePulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='PDB Keywords')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Keywords that will be in the header of the PDB file'

    self.pdbKeywordsMulti = MultiWidget(frameM, Entry, minRows=1, useImages=False, tipText=tipText,
                                        callback=None)
    self.pdbKeywordsMulti.grid(row=i, column=1, sticky='w')  
    self.pdbKeywordsMulti.bind('<Leave>', self.setPdbKeywords, '+')
    i += 1

    #label = Label(frameM, text='Structural Genomics\nProject Submission?')
    #label.grid(row=i, column=0, sticky='w')

    #self.sgCheck = CheckButton(frameM)
    #self.sgCheck.grid(row=i, column=1, sticky='w')
    #self.sgCheck.set(False)
    #i += 1

    """
    label = Label(frameM, text='Release Date - Coordinate Data')
    label.grid(row=i, column=0, sticky='w')

    coordDateTypes = ['release coordinates immediatedly',
                      'release coordinates upon publication (maximum hold period 1 year)',
                      'place coordinates on hold for one year']

    self.coordDatePulldown = PulldownList(frameM, texts=coordDateTypes,
                                          callback=self.setCoordDate)
    self.coordDatePulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Release Date - NMR Restraints')
    label.grid(row=i, column=0, sticky='w')

    restrDateTypes = ['release NMR restraints immediatedly (3 weeks after completion of deposition)',
                      'release NMR restraints data upon publication of paper',
                      'place one year hold on the NMR restraints']

    self.restrDatePulldown = PulldownList(frameM, texts=restrDateTypes,
                                          callback=self.setRestrDate)
    self.restrDatePulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Release Date - Other NMR Data')
    label.grid(row=i, column=0, sticky='w')

    nmrDateTypes = ['release other NMR data immediatedly (x weeks after completion of deposition)', # TODO - do we need to specify x?
                    'release other NMR data upon publication of paper',
                    'place one year hold on the other NMR data']

    self.nmrDatePulldown = PulldownList(frameM, texts=nmrDateTypes,
                                        callback=self.setNmrDate)
    self.nmrDatePulldown.grid(row=i, column=1, sticky='w')
    i += 1

    label = Label(frameM, text='Release Sequence Immediately?')
    label.grid(row=i, column=0, sticky='w')

    seqDateBools = ['yes', 'no']

    self.seqDatePulldown = PulldownList(frameM, texts=seqDateBools,
                                          callback=self.setSeqDate)
    self.seqDatePulldown.grid(row=i, column=1, sticky='w')
    i += 1
    """

    label = Label(frameM, text='Detailed Description')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Text for a detailed description of your study (like a mini abstract)'

    self.detailsText = Text(frameM, text='', width=32, height=5, tipText=tipText)
    self.detailsText.grid(row=i, column=1, sticky='ew')
    self.detailsText.bind('<Leave>',self.setDetails, '+')
    i += 1

    label = Label(frameM, text='Special Deposition\nInstructions')
    label.grid(row=i, column=0, sticky='w')

    tipText = 'Additional instructions or notes that you wish to convey to the curators about your deposition'

    self.processText = Text(frameM, text='', width=32, height=5, tipText=tipText)
    self.processText.grid(row=i, column=1, sticky='ew')  
    self.processText.bind('<Leave>', self.setProcessing, '+')
    i += 1

    div = LabelDivider(frameM, text='Completion Report', grid=(i,0), gridSpan=(1,2))
    i += 1

    headingList = ['CCPN Data Type','Comment']

    tipTexts = ['Indicates mandatory data items for PDB or BMRB depositions - clicking in any field will take you to the correct tab for editing the data', 'Has the mandatory data been completed? Finished data is shown in green. Clicking on any missing data (shown in red, orange, or yellow) will take you to the correct tab so that you can add this information']

    editWidgets = [None, None]
    editGetCallbacks = [self.remedyReport, self.remedyReport]
    editSetCallbacks = [None, None]

    self.reportTable = ScrolledMatrix(frameM, 
                                      multiSelect=False,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks, 
                                      editWidgets=editWidgets,
                                      headingList=headingList,
                                      tipTexts=tipTexts,
                                      callback=None, grid=(i,0), gridSpan=(1,2))   
    i += 1

    texts = ['Add New Deposition Entry','Remove Current Deposition Entry',
             'Export as NMR-STAR 3.1','Export PDB Coords'] #, 'Export CNS Distance Restraints', 'Make CCPN Project as tgz File']
    commands = [self.addNewEntry, self.deleteCurrentEntry,
                self.exportNmrStar31, self.exportPdb] #, self.exportCnsDistance, self.packageProject]

    tipTexts = ['Add a new "Entry" object that can be submitted to the PDB or BMRB',
                'Remove the currently selected "Entry"',
                'Export an NMR-STAR 3.1 file that can be uploaded in ADIT-NMR',
                'Export coordinates in PDB format for viewing with Rasmol, etc.'
                ]

    self.mainButtons = ButtonList(frameM, texts=texts, tipTexts=tipTexts, commands=commands)
    self.mainButtons.grid(row=i, column=0, columnspan=2, sticky='ew')
    self.mainButtons.buttons[0].config(bg=NICE_GREEN)

    #print 'COUNT: [%s]' % i

    i = 0

    #
    # People
    #

    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_columnconfigure(1, weight=1)
    frameA.grid_rowconfigure(1, weight=1)
    frameA.grid_rowconfigure(3, weight=1)

    div = LabelDivider(frameA, text='People')
    div.grid(row=0, column=0, columnspan=2, sticky='ew')

    self.familyNameEntry = Entry(self, returnCallback=self.setFamilyName)
    self.titlePulldown = PulldownList(self, texts=TITLES,
                                      callback=self.setPersonTitle)
    self.givenNameEntry = Entry(self, returnCallback=self.setGivenName)                                   
    self.initialsEntry = Entry(self, returnCallback=self.setInitials)
    self.famTitlePulldown = PulldownList(self, texts=FAM_TITLES,
                                         callback=self.setFamilyTitle)
    self.personGroupPulldown = PulldownList(self, callback=self.setPersonGroup)

    headingList = ['#','Entry Author?','Contact?','Primary\nContact?','Family Name','Title','Given Name',
                   'Initials','Family\nTitle (e.g. Jr.)','Group']

    tipTexts = ['Number of the person in list',
                'Is the person an author in the "Entry"? - toggle by clicking',
                'Is the person a contact person? - toggle by clicking',
                'Is the person the primary contact used by AutoDep? - toggle by clicking',
                'Family name of the person',
                'Professional title of the person such as Dr, Prof, Mr, etc.',
                'Given name of the person',
                "Any middle initials in the person's name",
                'Family title of the person like Sr., Jr., III, etc.',
                'Group name that the person is attached to - add new groups in the frame in the bottom left of this tab',
                ]

    editWidgets = [None, None, None, None, self.familyNameEntry, self.titlePulldown,
                    self.givenNameEntry, self.initialsEntry,
                    self.famTitlePulldown, self.personGroupPulldown]
    editGetCallbacks = [None, self.toggleAuthor, self.toggleContact, self.togglePrimary,
                        self.getFamilyName, self.getPersonTitle,
                        self.getGivenName, self.getInitials,
                        self.getFamilyTitle, self.getPersonGroup]
    editSetCallbacks = [None, None, None, None, self.setFamilyName, self.setPersonTitle,
                        self.setGivenName, self.setInitials,
                        self.setFamilyTitle, self.setPersonGroup]

    self.personTable = ScrolledMatrix(frameA, 
                                      multiSelect=False,
                                      headingList=headingList,
                                      tipTexts=tipTexts,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks, 
                                      editWidgets=editWidgets,
                                      callback=self.selectPersonTable)
    self.personTable.grid(row=1, column=0, columnspan=2, sticky='nsew')

    texts = ['Add Person','Remove Person','Set Person as Author',
             'Set Person As Contact', 'Set Person as Primary Contact',]

    tipTexts = ['Add a new person',
                'Remove the currently selected person',
                'Set the currently selected person as an "Entry" author',
                'Set the currently selected person as a contact person',
                'Set the currently selected person as the primary contact'
                ]

    commands = [self.addPerson, self.removePerson, self.setPersonAsAuthor,
                self.setPersonAsContact, self.setPersonAsPrimary]
    buttons = ButtonList(frameA, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, columnspan=2, sticky='ew')

    frame1 = LabelFrame(frameA, text='Group Information')
    frame1.grid(row=3, column=0, sticky='nsew')
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    #div = LabelDivider(frame1, text='Group Information', grid=(0,0))

    self.mailAddressEntry = Entry(self, returnCallback=self.setMailAddress)
    self.groupCityEntry = Entry(self, returnCallback=self.setGroupCity)
    self.stateProvEntry = Entry(self, returnCallback=self.setStateProv)
    self.countryEntry = Entry(self, returnCallback=self.setCountry)
    self.postalCodeEntry = Entry(self, returnCallback=self.setPostalCode)
    self.orgTypePulldown = PulldownList(self, texts=ORGN_TYPES,
                                        callback=self.setOrgType)

    headingList = ['Group Name','Mailing Address','City','State/Province',
                   'Postal\nCode', 'Country','Organisation\nType']

    tipTexts = ['Name of the group',
                'First few lines of the mailing address of the group',
                'City where the group is located',
                'State or Province such as Cambridgeshire, CA or ON',
                'Postal code or ZIP code',
                'Country where the group is located',
                'Type of organisation the group is attached to - from academic, government, commercial or other'
                ]

    editWidgets = [None, self.mailAddressEntry, self.groupCityEntry,
                   self.stateProvEntry, self.postalCodeEntry,
                   self.countryEntry, self.orgTypePulldown]

    editSetCallbacks = [None, self.setMailAddress, self.setGroupCity,
                        self.setStateProv, self.setPostalCode,
                        self.setCountry, self.setOrgType]

    editGetCallbacks = [None, self.getMailAddress, self.getGroupCity,
                        self.getStateProv, self.getPostalCode,
                        self.getCountry, self.getOrgType]

    self.groupTable = ScrolledMatrix(frame1,
                                     headingList=headingList,
                                     tipTexts=tipTexts,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks, 
                                     editWidgets=editWidgets,
                                     callback=self.selectGroupTable, grid=(0,0))

    texts = ['Add Group','Remove Group']

    tipTexts = ['Add a new group',
                'Remove the currently selected group']

    commands = [self.addGroup,  self.removeGroup]
    buttons = ButtonList(frame1, texts=texts, tipTexts=tipTexts, commands=commands, grid=(1,0))

    frame2 = LabelFrame(frameA, text='Contact Person Information')
    frame2.grid(row=3, column=1, sticky='nsew')
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)

    #div = LabelDivider(frame2, text='Contact Person Information', grid=(0,0))

    #self.emailAddressEntry = Entry(self, returnCallback=self.setEmailAddress)
    self.telephoneEntry = Entry(self, returnCallback=self.setTelephone)
    self.faxEntry = Entry(self, returnCallback=self.setFax)
    self.positionPulldown = PulldownList(self, texts=POS_ROLES,
                                         callback=self.setPosition)

    headingList = ['Family Name','Given Name','Email Address',
                   'Telephone','FAX','Position']

    tipTexts = ['Family name of the person',
                'Given name of the person',
                'Email address of the person - add with the button below',
                'Telephone number of the person',
                'FAX number of the person',
                'Role of the contact person in the project - principal investigator, responsible person, investigator'
                ]

    editWidgets = [None, None, None, #self.emailAddressEntry,
                   self.telephoneEntry, self.faxEntry, self.positionPulldown]

    editSetCallbacks = [None, None, None, #self.setEmailAddress,
                        self.setTelephone, self.setFax, self.setPosition]

    editGetCallbacks = [None, None, None, #self.getEmailAddress,
                        self.getTelephone, self.getFax, self.getPosition]

    self.addressTable = ScrolledMatrix(frame2,
                                       headingList=headingList,
                                       tipTexts=tipTexts,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks, 
                                       editWidgets=editWidgets,
                                       callback=self.selectAddressTable, grid=(0,0))

    #self.addressTable.doEditMarkExtraRules = self.doAddressTableEditMarkExtraRules

    texts = ['Add Address','Remove Address']

    tipTexts = ['Add a new email address for the person selected - the person does need to be set to being in a group using the two other frames in this tab',
                'Remove the email address and other information for the currently selected person - note it will also remove the link to their group']

    commands = [self.addAddress,  self.removeAddress]
    buttons = ButtonList(frame2, texts=texts, tipTexts=tipTexts, commands=commands, grid=(1,0))


    #
    # References
    #

    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(1, weight=1)
    frameB.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameB, text='Citations')
    div.grid(row=0, column=0, sticky='ew')

    self.citeDepositionPulldown = PulldownList(self, callback=self.setCiteDeposition)
    self.citationAuthorMulti = MultiWidget(self, PulldownList, minRows=1,
                                           callback=self.setCiteAuthors)
    self.citationEditorMulti = MultiWidget(self, PulldownList, minRows=1,
                                           callback=self.setCiteEditors)
    self.citeYearPulldown = PulldownList(self, callback=self.setCiteYear)
    self.citeTitleEntry = Text(self, text='', width=64, height=5)
    self.citeJournalPulldown = PulldownList(self, callback=self.setCiteJournal)
    self.citeIssueEntry = Entry(self, width=6, returnCallback=self.setCiteIssue)
    self.citeVolumeEntry = Entry(self, width=6, returnCallback=self.setCiteVolume)
    self.citeFirstPageEntry = IntEntry(self, width=6, returnCallback=self.setCiteFirstPage)
    self.citeLastPageEntry = IntEntry(self, width=6, returnCallback=self.setCiteLastPage)
    self.citeStatusPulldown = PulldownList(self, callback=self.setCiteStatus)
    self.citePubMedEntry = Entry(self, width=12, returnCallback=self.setCitePubMed)
    self.citeDoiEntry = Entry(self, width=12, returnCallback=self.setCiteDoi)
    self.citeKeywordsMulti = MultiWidget(self, Entry, minRows=2,
                                         callback=self.setCiteKeywords)
    self.citeDetailsEntry = Entry(self, width=12, returnCallback=self.setCiteDetails)

    self.bookTitleText = Text(self, text='', width=64, height=5)
    self.bookSeriesEntry = Entry(self, width=12, returnCallback=self.setBookSeries)
    self.bookPublisherEntry = Entry(self, width=12, returnCallback=self.setBookPublisher)
    self.citeCityEntry = Entry(self, width=12, returnCallback=self.setCiteCity)
    self.citeCountryEntry = Entry(self, width=12, returnCallback=self.setCiteCountry)
    self.thesisInstitutionEntry = Entry(self, width=12, returnCallback=self.setThesisInst)
    self.confTitleText = Text(self, text='', width=64, height=5)

    editWidgets = [None, None, self.citeDepositionPulldown,
                   self.citationAuthorMulti, self.citationEditorMulti,
                   self.citeYearPulldown, self.citeTitleEntry,
                   self.citeJournalPulldown, self.citeIssueEntry,
                   self.citeVolumeEntry, self.citeFirstPageEntry,
                   self.citeLastPageEntry, self.citeStatusPulldown,
                   self.citePubMedEntry, self.citeDoiEntry,
                   self.citeKeywordsMulti, self.citeDetailsEntry,
                   self.bookTitleText, self.bookSeriesEntry,
                   self.bookPublisherEntry, self.citeCityEntry,
                   self.citeCountryEntry, self.thesisInstitutionEntry,
                   self.confTitleText]

    editGetCallbacks = [None, None, self.getCiteDeposition,
                        self.getCiteAuthors, self.getCiteEditors,
                        self.getCiteYear, self.getCiteTitle,
                        self.getCiteJournal, self.getCiteIssue,
                        self.getCiteVolume, self.getCiteFirstPage,
                        self.getCiteLastPage, self.getCiteStatus,
                        self.getCitePubMed, self.getCiteDoi,
                        self.getCiteKeywords, self.getCiteDetails,
                        self.getBookTitle, self.getBookSeries,
                        self.getBookPublisher, self.getCiteCity,
                        self.getCiteCountry, self.getThesisInst,
                        self.getConfTitle]

    editSetCallbacks = [None, None, self.setCiteDeposition, 
                        self.setCiteAuthors, self.setCiteEditors,
                        self.setCiteYear, self.setCiteTitle,
                        self.setCiteJournal, self.setCiteIssue,
                        self.setCiteVolume, self.setCiteFirstPage,
                        self.setCiteLastPage, self.setCiteStatus,
                        self.setCitePubMed, self.setCiteDoi,
                        self.setCiteKeywords, self.setCiteDetails,
                        self.setBookTitle, self.setBookSeries,
                        self.setBookPublisher, self.setCiteCity,
                        self.setCiteCountry, self.setThesisInst,
                        self.setConfTitle]

    headingList = ['#','Type','Submit?',
                   'Authors','Editors',
                   'Year','Title',
                   'Journal\nAbbrev','Issue',
                   'Volume','First\npage',
                   'Last\nPage','Status',
                   'PubMed ID','DOI',
                   'Keywords', 'Details',
                   'Book\nTitle','Book\nSeries',
                   'Book\nPublisher','City','Country',
                   'Institution\ngranting Thesis',
                   'Conference\nTitle']

    tipTexts = ['Number of the citation in list',
                'Type of citation based on the button used to add it',
                'Is the citation to be submitted? Primary, other (for a secondary publication) or no',
                'Add authors for the citation based on people added in the preceding tab',
                'Add editors based on people added in the preceding tab',
                'Year of publication, only editable if the status of the citation is set to published',
                'Title of the citation (or title of any book chapters)',
                "The abbreviation of the journal's name (only for citations in journals)",
                'Issue number, only editable if the citation is from a published journal reference',
                'Volume number, only editable if the status of the citation is set to published',
                'First page number, only editable if the status of the citation is set to published',
                'Final page number, only editable if the status of the citation is set to published',
                'Publication status of the citation - from published, in press, in preparation',
                'PubMed number, only editable if the status of the citation is set to published',
                'Digital Object Identifier, only editable if the status of the citation is set to published',
                'Keywords used in the publication',
                'Additional information about the citation',
                'Title of the book (for book citations only)',
                'Book series (for book citations only)',
                'Publisher of the book (for book citations only)',
                'City of book publishers, institution of Thesis or conference location (not relevant for journal publications)',
                'Country of institution of Thesis or conference location (not needed for journals and books)',
                'Name of the institution granting a Thesis',
                'Name of conference for abstract citations'
                ]

    self.citationTable = ScrolledMatrix(frameB, 
                                        multiSelect=False,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        callback=self.selectCitation)
    self.citationTable.doEditMarkExtraRules = self.doCiteTableEditMarkExtraRules
    self.citationTable.grid(row=1, column=0, sticky='nsew')

    texts = ['Add Journal Paper','Add Book/Chapter','Add Thesis','Add Conference Paper','Remove Citation']

    tipTexts = ['Add a new journal citation',
                'Add a new book citation',
                'Add a new Thesis citation',
                'Add a new conference abstract citation',
                'Remove the currently selected citation'
                ]

    commands = [self.addJournal, self.addBook, self.addThesis, self.addConference, self.removeCitation]
    buttons = ButtonList(frameB, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, sticky='ew')

    div = LabelDivider(frameB, text='Software')
    div.grid(row=3, column=0, sticky='ew')

    self.softTasksMulti = MultiWidget(self, Entry, minRows=1,
                                      callback=self.setSoftTasks)
    #self.softVersionEntry = Entry(self, returnCallback=self.setSoftVersion)
    self.softVendorEntry = Entry(self, returnCallback=self.setSoftVendor)
    self.softVendorAddrEntry = Entry(self, returnCallback=self.setSoftVendorAddr)
    self.softVendorWebEntry = Entry(self, returnCallback=self.setSoftVendorWeb)
    self.softDetailsEntry = Entry(self, returnCallback=self.setSoftDetails)
    self.softMethodEntry = Entry(self, returnCallback=self.setSoftMethod)

    editWidgets = [None, None, self.softTasksMulti,
                   self.softVendorEntry, self.softVendorAddrEntry,
                   self.softVendorWebEntry, self.softDetailsEntry,
                   self.softMethodEntry]

    editGetCallbacks = [None, None, self.getSoftTasks,
                        self.getSoftVendor, self.getSoftVendorAddr,
                        self.getSoftVendorWeb, self.getSoftDetails,
                        self.getSoftMethod]

    editSetCallbacks = [None, None, self.setSoftTasks,
                        self.setSoftVendor, self.setSoftVendorAddr,
                        self.setSoftVendorWeb, self.setSoftDetails,
                        self.setSoftMethod]

    headingList = ['Name','Version','Tasks','Vendor/Author',
                   'Vendor\nAddress','Vendor\nWebsite','Details','Method']

    tipTexts = ['Name of the software used',
                'Version number of the software',
                'Choose a list of tasks that the software was used for during the study; such as peak picking, chemical shift assignment, data analysis, structure calculation, geometry optimisation, etc',
                'Authors or vendors who produced the software',
                'Mailing address for the software authors/vendors',
                'Web address for the software authors/vendors',
                'Additional information about the software',
                'Method that the software uses to perform the task; mainly for assignment/peak picking strategies, or for structure calculation methods such as simulated annealing, distance geometry, DGSA, molecular dynamics, matrix relaxation, torsion angle dynamics, etc.',
                ]

    self.softwareTable = ScrolledMatrix(frameB, 
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=True,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        callback=self.selectSoftware)
    self.softwareTable.grid(row=4, column=0, sticky='nsew')

    texts = ['Add Software','Remove Software']

    tipTexts = ['Add a new software and version of the software',
                'Remove the currently selected software'
                ]

    commands = [self.addSoftware, self.removeSoftware]
    buttons = ButtonList(frameB, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=5, column=0, sticky='ew')


    #
    # Nmr Data
    #

    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(1, weight=1)
    frameC.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameC, text='NMR Measurements')
    div.grid(row=0, column=0, sticky='ew')

    headingList = ['Data Type','Available','Selected','Notes']

    tipTexts = ['Type of NMR data such as shift lists, peak lists, relaxation data, etc.',
                'Number of lists in the CCPN projects for the specific NMR data type',
                'Number of lists selected for deposition to PDB/BMRB for the specific NMR data type',
                'Description of the NMR data type',
                ]

    editWidgets = []
    editGetCallbacks = []
    editSetCallbacks = []
    justification = ['center','center','center','left']

    self.nmrListsTable = ScrolledMatrix(frameC, 
                                        multiSelect=False,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        justifyList=justification, 
                                        callback=self.selectGenTable)

    self.nmrListsTable.grid(row=1, column=0, sticky='nsew')

    texts = ['Add All Available NMR data','Remove All NMR data']

    tipTexts = ['Rapidly select all available data lists for deposition to PDB/BMRB',
                'De-select all data lists for PDB/BMRB deposition - at least one chemical shift list is now mandatory for all depositions']

    commands = [self.addAllNmrLists, self.removeAllNmrLists]
    buttons = ButtonList(frameC, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, sticky='ew')

    self.submitText = 'Selected for Deposition'
    self.genDivider = LabelDivider(frameC, text='NMR Lists %s' % self.submitText)
    self.genDivider.grid(row=3, column=0, sticky='ew')

    self.nmrListNameEntry = Entry(self, text='', returnCallback=self.setNmrListName)
    self.nmrListSoftwarePulldown = PulldownList(self, callback=self.setNmrListSoftware)

    headingList = ['List','Name','Size','Submit?','Software']

    # Note: these are actually set in updateNmrListSelectTable -
    #   as headingList is also re-set there depending on the NMR list being selected.

    tipTexts = ['Indicate the type and name of each list, based on list type selectd in the top frame',
                'Name of list which can also be changed',
                'Number of data points in the list',
                'Is the list selected for deposition? - toggle by clicking',
                'Software used to make the list - the software will also need a method selected in the "References" tab'
                ]

    editWidgets = [None, self.nmrListNameEntry, None, None, self.nmrListSoftwarePulldown]
    editGetCallbacks = [None, self.getNmrListName, None, self.toggleNmrListSubmit, self.getNmrListSoftware]
    editSetCallbacks = [None, self.setNmrListName, None, None, self.setNmrListSoftware]
    self.nmrListSelectTable = ScrolledMatrix(frameC,
                                             multiSelect=True,
                                             headingList=headingList,
                                             tipTexts=tipTexts,
                                             editSetCallbacks=editSetCallbacks,
                                             editGetCallbacks=editGetCallbacks,
                                             editWidgets=editWidgets,
                                             callback=self.selectEditTable)
    self.nmrListSelectTable.doEditMarkExtraRules = self.doEditTableEditMarkExtraRules
    self.nmrListSelectTable.grid(row=4, column=0, sticky='nsew')

    texts = ['Select All','Select',
             'Unselect','Unselect All']

    tipTexts = ['Select all lists of the selected data type',
                'Submit the currently selected data list',
                'De-select all lists of the selected data type',
                'De-select the currently selected data list',
                ]

    commands = [self.submitAllNmrLists, self.submitSelectedNmrList,
                self.withdrawSelectedNmrList, self.withdrawAllNmrLists]
    self.editButtons = ButtonList(frameC, texts=texts, tipTexts=tipTexts, commands=commands)
    self.editButtons.grid(row=5, column=0, sticky='ew')


    #
    # Shifts
    #

    frameD.grid_columnconfigure(1, weight=1)
    frameD.grid_rowconfigure(1, weight=1)
    frameD.grid_rowconfigure(5, weight=1)

    # Shift references

    div = LabelDivider(frameD, text='Chemical Shift References')
    div.grid(row=0, column=0, columnspan=2, sticky='ew')

    self.isotopePulldown = PulldownList(self, texts=STANDARD_ISOTOPES,
                                        callback=self.setShiftRefIsotope)

    self.molNamePulldown = PulldownList(self, texts=SHIFT_REF_COMPOUNDS,
                                        callback=self.setShiftRefMolName)

    self.atomGroupEntry  = Entry(self, text='', width=8, returnCallback=self.setShiftRefAtomGroup)
    self.valueEntry      = FloatEntry(self, text=0.0, width=6, returnCallback=self.setShiftRefValue, formatPlaces=9)
    self.ratioEntry      = FloatEntry(self, text=0.0, width=6, returnCallback=self.setShiftRefRatio, formatPlaces=9)
    self.unitPulldown    = PulldownList(self,  texts=STANDARD_UNITS,
                                        callback=self.setShiftRefUnit)

    self.geometryEntry   = Entry(self, text='', returnCallback=self.setShiftRefGeometry)
    self.locationEntry   = Entry(self, text='', returnCallback=self.setShiftRefLocation)
    self.axisEntry       = Entry(self, text='', returnCallback=self.setShiftRefAxis)

    colHeadings      = ['#','Class','Isotope','Experiments',
                        'Mol. Name','Atom','Value','Unit',
                        'Ref Type','Indirect\nShift Ratio',
                        'Sample\nGeometry','Location','Axis']

    tipTexts = ['The serial number of the chemical shift reference specification',
                'Whether the chemical shift reference is internal or external to the sample',
                'The kind of nuclear isotope to which the reference applies'
                'The number of experiments in the project which use the shift reference specification',
                'The name of the molecule used to give a reference value to chemical shifts',
                'Which atom of the reference molecule provides the reference chemical shift value',
                'The reference value of the chemical shift for the specified atom',
                'Which measurement unit the chemical shift reference value is in; ppm, ppb or ppt',
                'Whether the chemical shift referencing is direct or indirect (and thus uses a shift ratio) - this field can be toggled',
                'The precise numeric ratio used to indirectly get the reference shift value of an isotope, given the direct measurement of a different isotope',
                'For external references, a description of the geometry of the container used to hold the reference compound - for example, cylindrical or spherical',
                'For external references, a description of the location of the reference',
                'For external references, orientation of the reference container with respect to external magnetic field - for example, parallel or perpendicular',
                ]

    editWidgets      = [None, None, self.isotopePulldown, None, self.molNamePulldown,
                        self.atomGroupEntry,self.valueEntry, self.unitPulldown,None,
                        self.ratioEntry,self.geometryEntry,self.locationEntry,self.axisEntry]

    editGetCallbacks = [None, None, self.getShiftRefIsotope, None, self.getShiftRefMolName,
                        self.getShiftRefAtomGroup,self.getShiftRefValue,self.getShiftRefUnit,self.toggleShiftRefType,
                        self.getShiftRefRatio ,self.getShiftRefGeometry,self.getShiftRefLocation,self.getShiftRefAxis]

    editSetCallbacks = [None, None, self.setShiftRefIsotope, None, self.setShiftRefMolName,
                        self.setShiftRefAtomGroup,self.setShiftRefValue,self.setShiftRefUnit,None,
                        self.setShiftRefRatio ,self.setShiftRefGeometry,self.setShiftRefLocation,self.setShiftRefAxis]

    self.shiftRefTable = ScrolledMatrix(frameD, multiSelect=True,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=colHeadings,
                                        tipTexts=tipTexts,
                                        callback=self.selectShiftRef)
    self.shiftRefTable.doEditMarkExtraRules = self.doShiftRefEditMarkExtraRules
    self.shiftRefTable.grid(row=1, column=0, columnspan=2, sticky='nsew')

    texts    = ['Add Standard\nIUPAC References','Add Internal\nReference','Add External\nReference','Remove\nSelected']

    tipTexts = ['Add standard IUPAC chemical shift references for 1H, 13C, 15N and 31P',
                'Add a new record for a chemical shift reference that is internal to the sample',
                'Add a new record for a chemical shift reference that is external to the sample',
                'Delete the selected chemical shift reference records',
                ]

    commands = [self.addStdIupacRef,self.addInternalShiftRef,self.addExternalShiftRef,self.removeShiftRefs]
    self.shiftRefButtons = ButtonList(frameD,texts=texts,tipTexts=tipTexts,expands=True,commands=commands)
    self.shiftRefButtons.grid(row=2, column=0, columnspan=2, sticky='ew')

    # Unusual shifts

    div = LabelDivider(frameD, text='Unusual Chemical Shifts')
    div.grid(row=3, column=0, columnspan=2, sticky='ew')

    label = Label(frameD, text='Shift List:')
    label.grid(row=4, column=0, sticky='w')

    tipText = 'Choose a chemical shift from the CCPN project'

    self.checkShiftListPulldown = PulldownList(frameD, tipText=tipText, callback=self.changeCheckShiftList)
    self.checkShiftListPulldown.grid(row=4, column=1, sticky='w')

    colHeadings = ['#','Isotope','Resonance','PPM',
                   'SD','BMRB\nMean','Random\nCoil']

    tipTexts = ['The serial number of the resonance',
                'The nuclear isotope type of the resonance',
                'The assignment annotation for the resonance',
                'Actual chemical shift of the resonance',
                'Standard deviation of shift in multiple experiments???', ##### is this right? #####
                'The average chemical shift for the atom type given in BioMagResBank data',
                'The sequence-adjusted random coil chemical shift value for the atom type',
                ]

    self.shiftCheckTable = ScrolledMatrix(frameD, multiSelect=False,
                                          headingList=colHeadings, tipTexts=tipTexts,
                                          callback=None)
    self.shiftCheckTable.grid(row=5, column=0, columnspan=2, sticky='nsew')


    #
    # Molecules
    #

    frameE.grid_columnconfigure(0, weight=1)
    frameE.grid_rowconfigure(1, weight=1)
    frameE.grid_rowconfigure(3, weight=1)

    div = LabelDivider(frameE, text='Molecular System Chains')
    div.grid(row=0, column=0, sticky='ew')

    headingList = ['Chain\nCode','Length','Molecule','Conformational\nIsomer',
                   'Chemical\nExchange State','Folded\nState']

    tipTexts = ['Code identifier for this chain',
                'Number of residues in chain',
                'Name of template molecule as shown in middle frame',
                'For biopolymers or ligands, does the indicated component of the system represent an observed conformational isomer of another component of the system. For example, component one may represent a polypeptide chain with proline 34 in the cis conformation and component two may represent the polypeptide chain with proline 34 in the trans conformation. In another example, component one may be alpha-D-glucose and component two beta-D-glucose',
                "A flag indicating whether the component of the system is in chemical exchange with another component of the system. Set as 'yes' if the component is in chemical exchange with another component in the system - for example ligands that are present in both the free and bound states",
                'An enumerated list of descriptive terms used to define the conformational state of the component of the assembly - for example, native, unfolded, etc.',
                ]

    self.chainConfIsomerPulldown    = PulldownList(self, callback=self.setChainConfIsomer)
    self.chainChemExchStatePulldown = PulldownList(self, callback=self.setChainChemExchState)
    self.chainPhysicalStatePulldown = PulldownList(self, callback=self.setChainPhysicalState)

    editWidgets = [None, None, None,
                   self.chainConfIsomerPulldown,
                   self.chainChemExchStatePulldown,
                   self.chainPhysicalStatePulldown]

    editGetCallbacks = [None, None, None,
                        self.getChainConfIsomer,
                        self.getChainChemExchState,
                        self.getChainPhysicalState,]

    editSetCallbacks = [None, None, None,
                        self.setChainConfIsomer,
                        self.setChainChemExchState,
                        self.setChainPhysicalState,]

    self.molSysTable = ScrolledMatrix(frameE, 
                                      multiSelect=False,
                                      headingList=headingList,
                                      tipTexts=tipTexts,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets,
                                      callback=self.selectChain)
    self.molSysTable.grid(row=1, column=0, sticky='nsew')

    div = LabelDivider(frameE, text='Molecules in Entry')
    div.grid(row=2, column=0, sticky='ew')

    headingList = ['Name', 'Molecule\nType', 'Number of\nResidues',
                   'Sequence', 'Is it\nParamagnetic?',
                   'Database\nName', 'Database\nAcc. Number',
                   'Database Details', 'Domain or\nFragment',
                   'EC Number', 'Mutations', 'Details']

    tipTexts = ['Name of template molecule',
                'Standard derived type of molecule - for example, protein, DNA, RNA, etc.',
                'Number of residues in this molecule',
                'One letter sequence of the molecule (shows the first 25 residues for longer polymers)',
                'Is this molecule paramagnetic? - derived boolean value',
                'Database name of cross reference (for example UniProt, GenBank, PDB, etc.), setting this will set a temporary value for the accession code, which should be set to the correct value',
                "Accession number for cross reference - it can only be set after adding a database name and will then be initialised to 'TmpAcc'",
                'Additional details pertaining to the database cross reference',
                'Give the name of the fragment if the protein or necleic acid is a part of a larger molecule with a well established biological function',
                'Enzyme Commission number if the molecule is an enzyme',
                'Specify any mutations that have been introduced into the molecular entity - for example, Y20A',
                'Additional details pertaining to this molecule or its sequence',
                ]

    self.moleculeDbNamePulldown = PulldownList(self, callback=self.setMoleculeDbName)
    self.moleculeDbRefEntry = Entry(self, returnCallback=self.setMoleculeDbRef)
    self.moleculeDbRefDetailsEntry = Entry(self, returnCallback=self.setMoleculeDbRefDetails)

    self.moleculeDomainEntry = Entry(self, returnCallback=self.setMoleculeDomain)
    self.moleculeECNumMulti= MultiWidget(self, Entry, minRows=1, callback=self.setMoleculeECnum)
    self.moleculeMutEntry = Entry(self, returnCallback=self.setMoleculeMutation)
    self.moleculeDetailsEntry = Entry(self, returnCallback=self.setMoleculeDetails)

    editWidgets = [None, None, None, None, None,
                   self.moleculeDbNamePulldown, self.moleculeDbRefEntry,
                   self.moleculeDbRefDetailsEntry,
                   self.moleculeDomainEntry, self.moleculeECNumMulti,
                   self.moleculeMutEntry, self.moleculeDetailsEntry]
    editGetCallbacks = [None, None, None, None, None,
                        self.getMoleculeDbName, self.getMoleculeDbRef,
                        self.getMoleculeDbRefDetails,
                        self.getMoleculeDomain, self.getMoleculeECnum,
                        self.getMoleculeMutation, self.getMoleculeDetails]
    editSetCallbacks = [None, None, None, None, None,
                        self.setMoleculeDbName, self.setMoleculeDbRef,
                        self.setMoleculeDbRefDetails,
                        self.setMoleculeDomain, self.setMoleculeECnum,
                        self.setMoleculeMutation, self.setMoleculeDetails]

    self.moleculeTable = ScrolledMatrix(frameE, 
                                        multiSelect=False,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        callback=self.selectMolecule)
    self.moleculeTable.grid(row=3, column=0, sticky='nsew')

    div = LabelDivider(frameE, text='Non-standard & Non-polymer Chem Comps')
    div.grid(row=4, column=0, sticky='ew')

    headingList = ['CcpCode', 'PDB Code', 'Mol. Type', 'Name']

    tipTexts = ['The CCPN code of the compound',
                'The PDB three-letter code of the compound',
                'Derived molecule type of chemical component - for example, protein, ???',
                'The long chemical name of the compound',
                ]

    self.unusualChemCompTable = ScrolledMatrix(frameE, 
                                               multiSelect=False,
                                               headingList=headingList,
                                               tipTexts=tipTexts,
                                               callback=None)
    self.unusualChemCompTable.grid(row=5, column=0, sticky='nsew')


    #
    # Biological sources
    #

    frameF.grid_columnconfigure(0, weight=1)
    frameF.grid_rowconfigure(1, weight=1)
    frameF.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameF, text='Experimental Sample Sources', grid=(0,0))

    self.expHostSciNamePulldown = PulldownList(self, callback=self.setExpHostSciName)
    self.expHostNameEntry = Entry(self, returnCallback=self.setExpHostName)
    self.expHostMolPulldown = PulldownList(self, callback=self.setExpHostMol)
    self.expHostProdMethodPulldown = PulldownList(self, callback=self.setExpHostProdMethod)
    self.expHostStrainEntry = Entry(self, returnCallback=self.setExpHostStrain)
    self.expHostVecTypePulldown = PulldownList(self, callback=self.setExpHostVecType)
    self.expHostVecNameEntry = Entry(self, returnCallback=self.setExpHostVecName)
    self.expHostCellLineEntry = Entry(self, returnCallback=self.setExpHostCellLine)
    self.expHostVariantEntry = Entry(self, returnCallback=self.setExpHostVariant)
    self.expHostDetailsEntry = Entry(self, returnCallback=self.setExpHostDetails)

    editWidgets = [self.expHostSciNamePulldown,self.expHostNameEntry,None,
                   self.expHostProdMethodPulldown,self.expHostStrainEntry,
                   self.expHostVecTypePulldown,self.expHostVecNameEntry,
                   self.expHostCellLineEntry,self.expHostVariantEntry,
                   self.expHostDetailsEntry]
    editGetCallbacks = [self.getExpHostSciName,self.getExpHostName,None,
                        self.getExpHostProdMethod,self.getExpHostStrain,
                        self.getExpHostVecType,self.getExpHostVecName,
                        self.getExpHostCellLine,self.getExpHostVariant,
                        self.getExpHostDetails]
    editSetCallbacks = [self.setExpHostSciName,self.setExpHostName,None,
                        self.setExpHostProdMethod,self.setExpHostStrain,
                        self.setExpHostVecType,self.setExpHostVecName,
                        self.setExpHostCellLine,self.setExpHostVariant,
                        self.setExpHostDetails]

    headingList = ['Scientific Name','Common Name','Molecule',
                   'Production\nMethod','Host strain','Vector\nType',
                   'Vector Name','Cell Line','Variant','Details']

    tipTexts = ["The standard scientific name for the source organism of this molecule - choose from a list of standard recombinant organisms or 'synthetic'",
                'The common name for the source organism of this molecule',
                'The name of the linked template polymeric molecule',
                'The method that was used to produce the organism - choose from a list including recombinant technology, chemical synthesis, etc.',
                'Strain information about the source organism',
                'The type of vector used such as plasmid, cosmid, etc.',
                'The specific name of the vector used if applicable',
                'The cell line that was used if applicable',
                'Variant details about the source organism if applicable',
                'Additional details about the production method or source organism',
                ]

    self.expSourceTable = ScrolledMatrix(frameF,
                                         headingList=headingList,
                                         tipTexts=tipTexts,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks, 
                                         editWidgets=editWidgets,
                                         callback=self.selectEntryMol, grid=(1,0))

    texts = ['Add Polymer Experimental Sources','Remove Experimental Source']

    tipTexts = ['Adds experimental and natural sources for all polymers in the molecular system',
                'Remove the currently selected experimental source']

    commands = [self.addExpSource, self.removeExpSource]
    buttons = ButtonList(frameF, texts=texts, tipTexts=tipTexts, commands=commands, grid=(2,0))

    div = LabelDivider(frameF, text='Biological & Sequence Origin', grid=(3,0))

    self.natSourceSciNamePulldown = PulldownList(self, callback=self.setNatSourceSciName)
    self.natSourceNameEntry = Entry(self, returnCallback=self.setNatSourceName)
    self.natSourceMolPulldown = PulldownList(self, callback=self.setNatSourceMol)
    self.natSourceOrgTypePulldown = PulldownList(self, callback=self.setNatSourceOrgType)
    self.natSourceStrainEntry = Entry(self, returnCallback=self.setNatSourceStrain)
    self.natSourceVariantEntry = Entry(self, returnCallback=self.setNatSourceVariant)
    self.natSourceGeneNameEntry = Entry(self, returnCallback=self.setNatSourceGeneName)
    self.natSourceCellLineEntry = Entry(self, returnCallback=self.setNatSourceCellLine)
    self.natSourceATCCNumEntry = Entry(self, returnCallback=self.setNatSourceATCCNum)
    self.natSourceOrganEntry = Entry(self, returnCallback=self.setNatSourceOrgan)
    self.natSourceTissueEntry = Entry(self, returnCallback=self.setNatSourceTissue)
    self.natSourceCellTypeEntry = Entry(self, returnCallback=self.setNatSourceCellType)
    self.natSourceOrganelleEntry = Entry(self, returnCallback=self.setNatSourceOrganelle)
    self.natSourceDetailsEntry = Entry(self, returnCallback=self.setNatSourceDetails)

    editWidgets = [self.natSourceSciNamePulldown,self.natSourceNameEntry,None,
                   self.natSourceOrgTypePulldown,self.natSourceStrainEntry,
                   self.natSourceVariantEntry,self.natSourceGeneNameEntry,
                   self.natSourceCellLineEntry,self.natSourceATCCNumEntry,
                   self.natSourceOrganEntry,self.natSourceTissueEntry,
                   self.natSourceCellTypeEntry,self.natSourceOrganelleEntry,
                   self.natSourceDetailsEntry]
    editGetCallbacks = [self.getNatSourceSciName,self.getNatSourceName,None,
                        self.getNatSourceOrgType,self.getNatSourceStrain,
                        self.getNatSourceVariant,self.getNatSourceGeneName,
                        self.getNatSourceCellLine,self.getNatSourceATCCNum,
                        self.getNatSourceOrgan,self.getNatSourceTissue,
                        self.getNatSourceCellType,self.getNatSourceOrganelle,
                        self.getNatSourceDetails]
    editSetCallbacks = [self.setNatSourceSciName,self.setNatSourceName,None,
                        self.setNatSourceOrgType,self.setNatSourceStrain,
                        self.setNatSourceVariant,self.setNatSourceGeneName,
                        self.setNatSourceCellLine,self.setNatSourceATCCNum,
                        self.setNatSourceOrgan,self.setNatSourceTissue,
                        self.setNatSourceCellType,self.setNatSourceOrganelle,
                        self.setNatSourceDetails]

    headingList = ['Scientific Name','Common Name','Molecule',
                   'Organism\nType','Strain','Variant','Gene Name',
                   'Cell Line','ATCC Number','Organ','Tissue',
                   'Cell Type','Organelle','Details']

    tipTexts = ['The standard scientific name for the biological or natural source of this molecule - choose from a list of the most common organisms in the PDB, which will also set the common name',
                'The common name for the biological or natural source organism of this molecule',
                'The name of the linked template polymeric molecule',
                'The type of organism that the molecule is from - for example, organism, virus, plasmid, no natural source, etc.',
                'Strain information about the source organism',
                'Variant details about the source organism if applicable',
                'Name of the gene that the protein is produced from',
                'The cell line that the protein is made from if applicable',
                'American Type Culture Collection number if sample comes from that organisation',
                'The organ in which the protein is located if applicable',
                'The tissue in which the protein is located if applicable',
                'The cell type that the protein is made from if applicable',
                'The organelle in which the protein is located if applicable',
                'Additional details about the host biological organism',
                ]

    self.natSourceTable = ScrolledMatrix(frameF,
                                         headingList=headingList,
                                         tipTexts=tipTexts,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks, 
                                         editWidgets=editWidgets,
                                         callback=self.selectNatSource, grid=(4,0))

    texts = ['Add Polymer Natural Sources','Remove Natural Source']

    tipTexts = ['Adds natural sources for all polymers in the molecular system (these are usually added when experimental sources are added)',
                'Remoevs the currently selected natural source']

    commands = [self.addNatSources, self.removeNatSource]
    buttons = ButtonList(frameF, texts=texts, tipTexts=tipTexts, commands=commands, grid=(5,0))


    #
    # Samples
    #

    frameG.grid_columnconfigure(0, weight=1)
    frameG.grid_rowconfigure(1, weight=1)
    frameG.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameG, text='Samples')
    div.grid(row=0, column=0, sticky='ew')

    #self.sampleStatePulldown = PulldownList(self, callback=self.setSampleState)

    self.sampleAmountEntry = FloatEntry(self, text='', width=12,
                                        returnCallback=self.setSampleAmount)

    self.sampleUnitPulldown = PulldownList(self, texts=SAMPLE_UNITS,
                                           callback=self.setSampleUnit)

    self.sampleIonicEntry = FloatEntry(self, text='', width=12,
                                       returnCallback=self.setSampleIonic)

    self.samplePhEntry = FloatEntry(self, text='', width=12,
                                       returnCallback=self.setSamplePh)

    self.sampleSolventPulldown = PulldownList(self, callback=self.setSampleSolvent)

    self.sampleDetailsEntry = Entry(self, text='', width=10,
                                 returnCallback=self.setSampleDetails)

    headingList = ['#','Name', #'Sample State',
                   #'Amount','Unit', 'Ionic\nStrength','pH',
                   'Solvent System','Details']

    tipTexts = ['The serial number of the sample',
                'The name used to label the sample',
                'The solvent composition of this sample - for example, 90% H2O/10% D2O, 100% D2O, etc.',
                'Additional details about the sample',
                ]

    editWidgets = [None, None, #self.sampleStatePulldown,
                   #self.sampleAmountEntry, self.sampleUnitPulldown,
                   #self.sampleIonicEntry, self.samplePhEntry,
                   self.sampleSolventPulldown, self.sampleDetailsEntry]

    editGetCallbacks = [None, None, #self.getSampleState,
                        #self.getSampleAmount, self.getSampleUnit,
                        #self.getSampleIonic, self.getSamplePh,
                        self.getSampleSolvent, self.getSampleDetails]

    editSetCallbacks = [None, None, #self.setSampleState,
                        #self.setSampleAmount, self.setSampleUnit,
                        #self.setSampleIonic, self.setSamplePh,
                        self.setSampleSolvent, self.setSampleDetails]

    self.sampleTable = ScrolledMatrix(frameG,
                                      multiSelect=False,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets,
                                      headingList=headingList,
                                      tipTexts=tipTexts,
                                      callback=self.selectSample)
    self.sampleTable.grid(row=1, column=0, sticky='nsew')

    texts = ['Add Sample','Remove Sample']

    tipTexts = ['Add a new sample to the CCPN project',
                'Remove the currently selected sample from the CCPN project']

    commands = [self.addSample, self.removeSample]
    buttons = ButtonList(frameG, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, sticky='ew')

    div = LabelDivider(frameG, text='Sample Components', grid=(3,0))

    #self.sampleCompMolPulldown = PulldownList(self, callback=self.setSampleCompMol)

    self.sampleCompIsotopePulldown = PulldownList(self, callback=self.setSampleCompIsotope)

    self.sampleCompConcEntry = FloatEntry(self, text='', width=12, formatPlaces=7,
                                          returnCallback=self.setSampleCompConc)

    self.sampleCompConcErrEntry = FloatEntry(self, text='', width=12, formatPlaces=7,
                                             returnCallback=self.setSampleCompConcErr)

    self.sampleCompConcUnitPulldown = PulldownList(self, texts=SAMPLE_COMP_UNITS,
                                           callback=self.setSampleUnit)

    headingList = ['Molecular Component','Isotopic Labelling','Concentration','Error','Unit']

    tipTexts = ['Common molecular name of the component from the sample selected in the top frame',
                'Isotopic labelling in standard IUPAC format - a common set is given in the pulldown menu or use the example provided to specify your own isotopic labelling',
                'Concentration value',
                'Error in concentration value',
                "Unit of concentration - from: 'kg/m3', 'M', 'm3/m3', 'mol/mol', 'kg/kg' - use 'm3/m3' to solute ratios",
                ]

    editWidgets = [None, self.sampleCompIsotopePulldown, self.sampleCompConcEntry,
                   self.sampleCompConcErrEntry, self.sampleCompConcUnitPulldown]

    editGetCallbacks = [None, self.getSampleCompIsotope, self.getSampleCompConc,
                        self.getSampleCompConcErr, self.getSampleCompConcUnit]

    editSetCallbacks = [None, self.setSampleCompIsotope, self.setSampleCompConc,
                        self.setSampleCompConcErr, self.setSampleCompConcUnit]

    self.sampleComponentTable = ScrolledMatrix(frameG,
                                               multiSelect=False,
                                               editSetCallbacks=editSetCallbacks,
                                               editGetCallbacks=editGetCallbacks,
                                               editWidgets=editWidgets,
                                               headingList=headingList,
                                               tipTexts=tipTexts,
                                               callback=self.selectSampleComponent, grid=(4,0))

    texts = ['Add Polymer Sample Components','Add Other Sample Component','Remove Sample Component']

    tipTexts = ['Add the names of all polymeric components from your molecular system',
                'Add the names of other components in this sample',
                'Remove the currently selected component',
                ]

    commands = [self.addPolySampleComponents, self.addOtherSampleComponent, self.removeSampleComponent]
    buttons = ButtonList(frameG, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=5, column=0, sticky='ew')


    #
    # NMR Conditions
    #

    frameH.grid_columnconfigure(0, weight=1)
    frameH.grid_rowconfigure(1, weight=1)
    frameH.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameH, text='NMR Condition Sets')
    div.grid(row=0, column=0, sticky='ew')

    self.sampleConditionSetEntry = Entry(self, text='', width=10,
                                         returnCallback=self.setSampleConditionSet)

    self.sampleConditionSetDetailsEntry = Entry(self, text='', width=10,
                                                returnCallback=self.setSampleConditionSetDetails)

    headingList = ['#','Condition Set','Details']

    tipTexts = ['The serial number of the condition set',
                'The name used to label the condition set',
                'Additional details about the condition set',
                ]

    editWidgets = [None,self.sampleConditionSetEntry,self.sampleConditionSetDetailsEntry]

    editGetCallbacks = [None,self.getSampleConditionSet,self.getSampleConditionSetDetails]

    editSetCallbacks = [None,self.setSampleConditionSet,self.setSampleConditionSetDetails]

    self.nmrCondSetTable = ScrolledMatrix(frameH,
                                          headingList=headingList,
                                          tipTexts=tipTexts,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          editWidgets=editWidgets,
                                          callback=self.selectSampleConditionSet)

    self.nmrCondSetTable.grid(row=1, column=0, sticky='nsew')

    texts = ['Add Condition Set','Remove Condition Set']

    tipTexts = ['Add a new condition set to the CCPN project',
                'Remove the currently select condition set from the CCPN project',
                ]

    commands = [self.addSampleConditionSet, self.removeSampleConditionSet]

    buttons = ButtonList(frameH, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, sticky='ew')

    div = LabelDivider(frameH, text='NMR Conditions')
    div.grid(row=3, column=0, sticky='ew')

    self.sampleConditionPulldown = PulldownList(self, callback=self.setSampleCondition)

    self.sampleConditionValueEntry = FloatEntry(self, text='', width=12,
                                                returnCallback=self.setSampleConditionValue)

    self.sampleConditionErrorEntry = FloatEntry(self, text='', width=12,
                                                returnCallback=self.setSampleConditionError)

    self.sampleConditionUnitEntry = Entry(self, text='', width=10,
                                          returnCallback=self.setSampleConditionUnit)

    headingList = ['Condition','Value','Error','Unit']

    tipTexts = ['The name of the variable used to define a specific sample condition - for example, temperature, pH, etc.',
                'The experimental value of the condition selected',
                'The experiemental error in the condition selected',
                'The unit used to measure the condition - for example, K for temperature'
                ]

    editWidgets      = [None,self.sampleConditionValueEntry,
                        self.sampleConditionErrorEntry,self.sampleConditionUnitEntry]

    editGetCallbacks = [None,self.getSampleConditionValue,
                        self.getSampleConditionError,self.getSampleConditionUnit]

    editSetCallbacks = [None,self.setSampleConditionValue,
                        self.setSampleConditionError,self.setSampleConditionUnit]

    self.nmrCondTable = ScrolledMatrix(frameH,
                                       headingList=headingList,
                                       tipTexts=tipTexts,
                                       editSetCallbacks=editSetCallbacks,
                                       editGetCallbacks=editGetCallbacks, 
                                       editWidgets=editWidgets,
                                       callback=self.selectSampleCondition)

    self.nmrCondTable.grid(row=4, column=0, sticky='nsew')

    texts = ['Add Standard Conditions','Add Condition','Remove Condition']

    tipTexts = ['This button adds the four standard condtions: temperature, pH, pressure and ionic strength',
                'Add other types of condition to the set',
                'Remove currently selected condition from the set',
                ]
                

    commands = [self.addStdSampleCondition, self.addSampleCondition, self.removeSampleCondition]

    buttons = ButtonList(frameH, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=5, column=0, sticky='ew')


    #
    # Instruments
    #

    frameI.grid_columnconfigure(0, weight=1)
    frameI.grid_rowconfigure(1, weight=1)
    frameI.grid_rowconfigure(4, weight=1)

    div = LabelDivider(frameI, text='NMR Spectrometers')
    div.grid(row=0, column=0, sticky='ew')

    self.spectrometerNameEntry = Entry(self, text='', width=10,
                                       returnCallback=self.setSpectrometerName)

    self.spectrometerFreqEntry = FloatEntry(self, text='', width=12,
                                            returnCallback=self.setSpectrometerFreq)

    self.spectrometerManufacturerPulldown = PulldownList(self, callback=self.setSpectrometerManufacturer)

    self.spectrometerModelPulldown = PulldownList(self, callback=self.setSpectrometerModel)

    self.spectrometerSerialEntry = Entry(self, text='', width=10,
                                         returnCallback=self.setSpectrometerSerial)

    self.spectrometerDetailsEntry = Entry(self, text='', width=10,
                                          returnCallback=self.setSpectrometerDetails)

    headingList = ['#','Name','Nominal Freq.','Proton Freq. (MHz)',
                   'Manufacturer','Model','Serial #','Details']

    tipTexts = ['Serial number of the NMR spectrometer specification',
                'A name for the spectrometer, for graphical displays',
                "The rounded spectrometer frequency, from the 1H resonance frequency in MHz, used in textual description - for example. '500', '900'",
                "The actual numeric magnetic field strength expressed as a 1H resonance frequency in MHz (also sets and rounds the nominal value) - for example '500.013'",
                'Manufacturer name - choose from a list including Bruker, Varian, GEOL, Custom, etc.',
                "The manufacturer's definition of the spectrometer model - choose from a list based on the manufacturer",
                "The manufacturer's serial number for the specific NMR spectrometer",
                'A user-specified textual comment about the NMR spectrometer',
                ]

    editWidgets      = [None,self.spectrometerNameEntry,
                        None,self.spectrometerFreqEntry,
                        self.spectrometerManufacturerPulldown,
                        self.spectrometerModelPulldown,
                        self.spectrometerSerialEntry,
                        self.spectrometerDetailsEntry]

    editGetCallbacks = [None,self.getSpectrometerName,
                        None,self.getSpectrometerFreq,
                        self.getSpectrometerManufacturer,
                        self.getSpectrometerModel,
                        self.getSpectrometerSerial,
                        self.getSpectrometerDetails]

    editSetCallbacks = [None,self.setSpectrometerName,
                        None,self.setSpectrometerFreq,
                        self.setSpectrometerManufacturer,
                        self.setSpectrometerModel,
                        self.setSpectrometerSerial,
                        self.setSpectrometerDetails]

    self.spectrometerMatrix = ScrolledMatrix(frameI,
                                             editSetCallbacks=editSetCallbacks,
                                             editGetCallbacks=editGetCallbacks,
                                             editWidgets=editWidgets,
                                             headingList=headingList,
                                             tipTexts=tipTexts,
                                             callback=self.selectSpectrometer)

    self.spectrometerMatrix.grid(row=1, column=0, sticky='nsew')

    texts    = ['Add Spectrometer','Remove Spectrometer']

    tipTexts = ['Add a new NMR spectrometer specification to the CCPN project',
                'Delete the selected NMR spectrometer specification from the CCPN project']

    commands = [self.addSpectrometer,self.removeSpectrometer]
    self.spectrometerButtons = ButtonList(frameI,texts=texts, tipTexts=tipTexts,expands=True, commands=commands)
    self.spectrometerButtons.grid(row=2, column=0, sticky='ew')

    div = LabelDivider(frameI, text='NMR Probes')
    div.grid(row=3, column=0, sticky='ew')

    self.probeNameEntry = Entry(self, text='', width=10,
                                returnCallback=self.setProbeName)

    self.probeTypePulldown = PulldownList(self, texts=NMR_PROBE_TYPES,
                                          callback=self.setProbeType)
    
    self.probeManufacturerEntry = Entry(self, text='', width=10,
                                        returnCallback=self.setProbeManufacturer)

    self.probeModelEntry = Entry(self, text='', width=10,
                                 returnCallback=self.setProbeModel)

    self.probeSerialEntry = Entry(self, text='', width=10,
                                  returnCallback=self.setProbeSerial)
                                       
    self.probeDiameterEntry = FloatEntry(self, text=0.0, width=10,
                                         returnCallback=self.setProbeDiameter)
                                   
    self.probeDetailsEntry = Entry(self, text='', width=10,
                                   returnCallback=self.setProbeDetails)

    headingList = ['#','Name','Type','Manufacturer','Model','Serial #','Diameter (cm)','Details']

    tipTexts = ['Serial number of NMR probe specification',
                'The name of the probe for graphical representation',
                'A classification for the kind of probe used - for example, liquid, solid, nano, flow or MAS',
                'Manufacturer name',
                "The manufacturer's definition of the probe model",
                "The manufacturer's serial number for the specific NMR probe",
                'The probe diameter in cm',
                'A user-specified textual comment about the probe',
                ]

    editWidgets      = [None,self.probeNameEntry, self.probeTypePulldown,
                        self.probeManufacturerEntry, self.probeModelEntry, self.probeSerialEntry,
                        self.probeDiameterEntry, self.probeDetailsEntry]
                        
    editGetCallbacks = [None,self.getProbeName,self.getProbeType,
                        self.getProbeManufacturer,self.getProbeModel,self.getProbeSerial,
                        self.getProbeDiameter,self.getProbeDetails]
    
    editSetCallbacks = [None,self.setProbeName,self.setProbeType,
                        self.setProbeManufacturer,self.setProbeModel,self.setProbeSerial,
                        self.setProbeDiameter,self.setProbeDetails]
                        
    self.probeMatrix = ScrolledMatrix(frameI,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets,
                                      headingList=headingList,
                                      tipTexts=tipTexts,
                                      callback=self.selectProbe)
                                      
    self.probeMatrix.grid(row=4, column=0, sticky='nsew')

    texts    = ['Add Probe','Remove Probe']

    tipTexts = ['Add a new NMR probe specification to the CCPN project',
                'Remove the selected NMR probe specification from the CCPN project']

    commands = [self.addProbe,self.removeProbe]
    self.probeButtons = ButtonList(frameI,texts=texts, tipTexts=tipTexts, expands=True, commands=commands)
    self.probeButtons.grid(row=5, column=0, sticky='ew')


    # Experiments

    frameJ.grid_columnconfigure(0, weight=1)
    frameJ.grid_rowconfigure(1, weight=1)

    div = LabelDivider(frameJ, text='Experiments')
    div.grid(row=0, column=0, sticky='ew')

    self.experimentNameEntry = Entry(self, text='', width=10,
                                     returnCallback=self.setExperimentName)

    self.experimentStatePulldown = PulldownList(self, callback=self.setExperimentState)

    self.experimentSamplePulldown = PulldownList(self, callback=self.setExperimentSample)

    self.experimentSampleCondPulldown = PulldownList(self, callback=self.setExperimentSampleCond)

    self.experimentSpecPulldown = PulldownList(self, callback=self.setExperimentSpectrometer)

    self.experimentProbePulldown = PulldownList(self, callback=self.setExperimentProbe)

    self.expShiftRefSelect = MultiWidget(self, CheckButton, callback=self.setExpShiftRefs,
                                           minRows=0, useImages=False)

    self.expChemShiftListPulldown = PulldownList(self, callback=self.setExpChemShiftList)

    headingList = ['Name','Submit?','Systematic\nName','Sample\nState','Sample','Condition\nSet','Chemical\nShift Ref.',
                   'Raw\ndata?','Spectrometer','Probe','Chemical\nShift List']

    tipTexts = ['The textual name for the experiment, for graphical displays',
                'Is this experiment to be submitted? - toggle by clicking (note if the experiment is part of a selected peak list, it will automatically be submitted)',
                'The name defined by the CCPN data model for this experiment or the user-defined name if no reference NMR experiment has been defined',
                'The state that best describes the sample and any molecular ordering; liquid (solution), solid, powder, ordered or crystalline',
                'Choose which sample was used to perform the experiment from the list added to the CCPN project in the "Samples" tab',
                'Choose which sample condition set was used to perform the experiment from the list added to the CCPN project in the "NMR Conditions" tab',
                'Choose multiple chemical shift references were used for any chemical shift list that the experiment was used to compile from the list added to the CCPN project in the "Shifts" tab',
                'Is the raw data available?',
                'Choose which NMR spectrometer was used to perform the experiment from the list added to the CCPN project in the "Instruments" tab',
                'Choose which NMR probe was used to perform the experiment from the list added to the CCPN project in the "Instruments" tab',
                'Choose a chemical shift list that the NMR experiment was used to compile, if applicable - the shift list should be selected for submission in the "NMR Data" tab',
                ]


    editWidgets = [self.experimentNameEntry, None, None,
                   self.experimentStatePulldown, self.experimentSamplePulldown,
                   self.experimentSampleCondPulldown, self.expShiftRefSelect,
                   None, self.experimentSpecPulldown, self.experimentProbePulldown,
                   self.expChemShiftListPulldown]
    editGetCallbacks = [self.getExperimentName, self.toggleExperiment, None,
                        self.getExperimentState, self.getExperimentSample,
                        self.getExperimentSampleCond, self.getExpShiftRefs,
                        None, self.getExperimentSpectrometer, self.getExperimentProbe,
                        self.getExpChemShiftList]
    editSetCallbacks = [self.setExperimentName, None, None,
                        self.setExperimentState, self.setExperimentSample,
                        self.setExperimentSampleCond, self.setExpShiftRefs,
                        None, self.setExperimentSpectrometer, self.setExperimentProbe,
                        self.setExpChemShiftList]

    self.experimentTable = ScrolledMatrix(frameJ,
                                          multiSelect=False,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          editWidgets=editWidgets,
                                          headingList=headingList,
                                          tipTexts=tipTexts,
                                          callback=self.selectExperiment)
    self.experimentTable.grid(row=1, column=0, sticky='nsew')

    texts = ['Add Experiment','Remove Experiment','Submit All','Submit Minimal']

    tipTexts = ['Add a new NMR experiment record to the CCPN project',
                'Remove the currently selected NMR experiment from the CCPN project',
                'Add all NMR experiments in the list to the "Entry" deposition object',
                'Just submit a minimal list of NMR experiments - any experiments used in submitted peak lists will remain to be submitted',
                ]

    commands = [self.addExperiment, self.removeExperiment,
                self.submitAllExperiments, self.submitMinExperiments]
    buttons = ButtonList(frameJ, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, sticky='ew')


    #
    # Structures
    #

    frameK.grid_columnconfigure(1, weight=1)
    frameK.grid_rowconfigure(1, weight=1)
    frameK.grid_rowconfigure(4, weight=1)
    frameK.grid_rowconfigure(7, weight=1)

    div = LabelDivider(frameK, text='Structure Generations')
    div.grid(row=0, column=0, columnspan=2, sticky='ew')

    self.strucGenNameEntry = Entry(self, returnCallback=self.setStrucGenName)
    self.strucGenTypePulldown = PulldownList(self, callback=self.setStrucGenType)
    self.strucGenDetailsEntry = Entry(self, returnCallback=self.setStrucGenDetails)
    self.strucGenSoftwarePulldown = PulldownList(self, callback=self.setStrucGenSoftware)

    self.strucEnsPulldown = PulldownList(self, callback=self.setStrucEns)
    self.strucGenConstSetPulldown = PulldownList(self, callback=self.setStrucGenConstSet)

    editWidgets = [None, self.strucGenNameEntry, None,
                   self.strucGenTypePulldown, self.strucGenDetailsEntry,
                   self.strucGenSoftwarePulldown, self.strucEnsPulldown,
                   self.strucGenConstSetPulldown]

    editGetCallbacks = [None, self.getStrucGenName, self.toggleStrucGen,
                        self.getStrucGenType, self.getStrucGenDetails,
                        self.getStrucGenSoftware, self.getStrucEns,
                        self.getStrucGenConstSet]

    editSetCallbacks = [None, self.setStrucGenName, None,
                        self.setStrucGenType, self.setStrucGenDetails,
                        self.setStrucGenSoftware, self.setStrucEns,
                        self.setStrucGenConstSet]

    headingList = ['Structure Generation #', 'Name', 'Submit?', 'Generation Type',
                   'Details', 'Software', 'Ensemble #', 'Constraint Set #']

    tipTexts = ['The serial number of the structure generation',
                'The name used to label the structure generation',
                'Is this structure generation to be submitted? - toggle by clicking',
                'Type of structure generation - denovo or refinement',
                'Additional details about the structure generation',
                'Software used to make the structure ensemble in this structure generation - the software will also need a method selected in the "References" tab',
                'Select from the structure ensembles in the CCPN project (shown in the middle frame) that was calculated using this structure generation',
                'Select from the constraint sets in the CCPN project (shown in the bottom frame) that was used to calculate the structure ensemble in this structure generation',
                ]

    self.strucGenTable = ScrolledMatrix(frameK, 
                                        multiSelect=False,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        callback=self.selectStrucGen)
    self.strucGenTable.grid(row=1, column=0, columnspan=2, sticky='nsew')

    texts = ['Add Structure Generation','Remove Structure Generation']

    tipTexts = ['Add a new structure generation to the CCPN project',
                'Remove the currently selected structure generation from the CCPN project',
                ]

    commands = [self.addStrucGen, self.removeStrucGen]
    buttons = ButtonList(frameK, texts=texts, tipTexts=tipTexts, commands=commands)
    buttons.grid(row=2, column=0, columnspan=2, sticky='ew')

    div = LabelDivider(frameK, text='Structure Ensembles')
    div.grid(row=3, column=0, columnspan=2, sticky='ew')

    self.strucEnsDetailsEntry = Entry(self, returnCallback=self.setStrucEnsDetails)
    self.strucEnsCalcEntry = Entry(self, returnCallback=self.setStrucEnsCalc)
    self.strucEnsCriteriaEntry = Entry(self, returnCallback=self.setStrucEnsCriteria)
    self.strucEnsReprEntry = Entry(self, returnCallback=self.setStrucEnsRepr)
    self.strucEnsReprCriEntry = Entry(self, returnCallback=self.setStrucEnsReprCri)

    editWidgets = [None, self.strucEnsDetailsEntry, None, None, None, self.strucEnsCalcEntry,
                   self.strucEnsCriteriaEntry, self.strucEnsReprEntry, self.strucEnsReprCriEntry]
    
    editGetCallbacks = [None, self.getStrucEnsDetails, None, None, None, self.getStrucEnsCalc,
                        self.getStrucEnsCriteria, self.getStrucEnsRepr,
                        self.getStrucEnsReprCri]

    editSetCallbacks = [None, self.setStrucEnsDetails, None, None, None, self.setStrucEnsCalc,
                        self.setStrucEnsCriteria, self.setStrucEnsRepr,
                        self.setStrucEnsReprCri]

    headingList = ['Structure Ensemble #', 'Details', 'Models', 'Chains', 'Residues',
                   'Calculated\nModels', 'Selection\nCriteria',
                   'Representative\nModel', 'Selection Criteria\nfor Representative']

    tipTexts = ['The serial number of the structure ensemble',
                'Additional details about this structure ensemble',
                'Number of models selected in this structure ensemble',
                'A list of chains in the molecular system that are represented in the structure ensemble',
                'The number of residues summed over all chains represented in the structure ensemble',
                'Number of models calculated in this structure ensemble',
                'Criteria used to select models to represent the structure ensemble from all calculated models',
                'The representative model for the structure ensemble',
                'Criteria used to select the representative model in the structure ensemble',
                ]

    self.strucEnsTable = ScrolledMatrix(frameK, 
                                        multiSelect=False,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        tipTexts=tipTexts,
                                        callback=self.selectStrucEns)
    self.strucEnsTable.grid(row=4, column=0, columnspan=2, sticky='nsew')

    div = LabelDivider(frameK, text='Constraint Sets')
    div.grid(row=5, column=0, columnspan=2, sticky='ew')

    label = Label(frameK, text='Constraint Set:')
    label.grid(row=6, column=0, sticky='w')

    tipText = 'Choose a constraint set from the CCPN project (there may be multiple constraint lists in each set)'

    self.constraintSetPulldown = PulldownList(frameK, tipText=tipText, callback=self.changeConstraintSet)
    self.constraintSetPulldown.grid(row=6, column=1, sticky='w')

    self.constListNameEntry = Entry(self, returnCallback=self.setConstListName)
    self.constListDetailsEntry = Entry(self, returnCallback=self.setConstListDetails)

    editWidgets = [None, None, self.constListNameEntry, self.constListDetailsEntry, None]

    editGetCallbacks = [None, None, self.getConstListName, self.getConstListDetails, None]

    editSetCallbacks = [None, None, self.setConstListName, self.setConstListDetails, None]

    headingList = ['Constraint List #', 'Constraint Type', 'Name', 'Details', 'Number of Constraints']

    tipTexts = ['The serial number of the constraint list',
                'The type of constraints in this list - for example, distances, NOEs, dihedral angles, RDCs, etc',
                'The name used to label this constraint list',
                'Additional details about this constraint list',
                'The number of constraints found in this list',
                ]

    self.constraintListsTable = ScrolledMatrix(frameK, 
                                               multiSelect=False,
                                               editSetCallbacks=editSetCallbacks,
                                               editGetCallbacks=editGetCallbacks,
                                               editWidgets=editWidgets,
                                               headingList=headingList,
                                               tipTexts=tipTexts,
                                               callback=self.selectConstraintList)
    self.constraintListsTable.grid(row=7, column=0, columnspan=2, sticky='nsew')


    #
    # Back to main
    #

    self.administerNotifiers(self.registerNotify)


  # Overall update functions

  def administerNotifiers(self, notifyFunc):

    for func in ('','__init__','delete'):
      for clazz in ('ccp.nmr.NmrEntry.Entry',):
        notifyFunc(self.updateEntryAfter,clazz, func)

    # People

    for func in ('','__init__','delete','setPersonInGroup',
                 'setFamilyName', 'setTitle', 'setGivenName',
                 'setMiddleInitials', 'setFamilyTitle'):
      notifyFunc(self.updatePersonAfter, 'ccp.general.Affiliation.Person', func)

    notifyFunc(self.updatePersonAfter, 'ccp.general.Affiliation.PersonInGroup', 'setGroup')

    # Groups

    for func in ('','__init__','delete','setName',):
      notifyFunc(self.updateGroupAfter, 'ccp.general.Affiliation.Group', func)

    for func in ('','__init__','delete','setName','setAddresses','setCity',
                 'setCountry','setPostalCode','setOrganisationType'):
      notifyFunc(self.updateGroupAfter, 'ccp.general.Affiliation.Organisation', func)

    # Contact people

    for func in ('','__init__','delete','setGroup','setEmailAddress',
                 'setPhoneNumbers','setFaxNumber','setPosition'):
      notifyFunc(self.updateAddressAfter, 'ccp.general.Affiliation.PersonInGroup', func)

    # Citations

    for func in ('','__init__','delete'):
      for clazz in ('ccp.general.Citation.JournalCitation','ccp.general.Citation.BookCitation',
                    'ccp.general.Citation.ThesisCitation','ccp.general.Citation.ConferenceCitation'):
        notifyFunc(self.updateCitations, clazz,func)

    # Software

    for func in ('','__init__','delete', 'setVersion', 'setTasks', 'setDetails',
                 'setVendorName', 'setVendorAddress', 'setVendorWebAddress', 'setMethods'):
      notifyFunc(self.updateSoftwareAfter, 'ccp.general.Method.Software',func)

    # Molecules & chains
    
    notifyFunc(self.updateMolDataAfter,'ccp.nmr.NmrEntry.Entry', 'setMolSystem')

    for func in ('','__init__','delete',
                 'setConformationalIsomer',
                 'setPhysicalState',
                 'setChemExchangeState',):
      notifyFunc(self.updateMolDataAfter,'ccp.molecule.MolSystem.Chain', func)

    for func in ('','__init__','delete','newAlignment','setDetails','addApplicationData'):
      notifyFunc(self.updateMolDataAfter,'ccp.molecule.Molecule.Molecule', func)

    for func in ('','__init__','delete','newName'):
      notifyFunc(self.updateMolDataAfter,'ccp.molecule.Molecule.MoleculeSysName', func)

    # Cross refs

    for func in ('','__init__','delete','newDatabase','findFirstDatabase'):
      notifyFunc(self.updateMolDataAfter,'memops.Implementation.MemopsRoot', func)

    for func in ('','__init__','delete','newEntry'):
      notifyFunc(self.updateMolDataAfter,'ccp.general.DbRef.Database', func)

    for func in ('','__init__','delete','setName','setCode','setDetails'):
      for clazz in ('ccp.general.DbRef.Database','ccp.general.DbRef.Entry'):
        notifyFunc(self.updateMolDataAfter,clazz, func)

    # Shift references

    for func in ('','__init__','delete','setAtomGroup','setIndirectShiftRatio',
                 'setValue','setIsotopeCode','setMolName','setReferenceType',
                 'setUnit','setExperiments','addExperiment','removeExperiment'):
      for clazz in ('ccp.nmr.Nmr.ExternalShiftReference','ccp.nmr.Nmr.InternalShiftReference'):
        notifyFunc(self.updateShiftRefs,clazz, func)

    for func in ('','__init__','delete','setShiftReferences',
                 'addShiftReference','removeShiftReference'):
      notifyFunc(self.updateShiftRefs,'ccp.nmr.Nmr.Experiment', func)

    # Sources

    for func in ('','__init__','delete','setMolecule','setProductionMethod','setSourceType','setVectorType'):
      notifyFunc(self.updateSourcesAfter, 'ccp.nmr.NmrEntry.EntryMolecule', func)

    for func in ('','__init__','delete','setScientificName','setOrganismName','setStrain',
                 'setPlasmid','setVariant','setGeneMnemonic','setDetails','setCellLine',
                 'setAtccNumber','setCellType','setOrgan','setOrganelle','setTissue'):
      notifyFunc(self.updateSourcesAfter,'ccp.general.Taxonomy.NaturalSource', func)

    # Samples

    for func in ('','__init__','delete','setName',
                 'setInitialAmount', 'setAmountUnit',
                 'setIonicStrength', 'setPh', 'setDetails'):
      notifyFunc(self.updateSamplesAfter,'ccp.lims.Sample.Sample', func)

    # Sample components

    for func in ('','__init__','delete','setName','setDetails'):
      notifyFunc(self.updateSampleComponentsAfter,'ccp.lims.RefSampleComponent.MolComponent', func)

    for func in ('','__init__','delete','setRefComponent','setConcentration',
                 'setConcentrationError','setConcentrationUnit'):
      notifyFunc(self.updateSampleComponentsAfter,'ccp.lims.Sample.SampleComponent', func)

    #for func in ('','__init__','delete','newResLabelFraction', 'newUniformAtomLabel', 'newSingleAtomLabel'):
    #  notifyFunc(self.updateSampleComponentsAfter,'ccp.molecule.LabeledMolecule.ResLabel', func)

    #for func in ('','__init__','delete','setWeight'):
    #  for clazz in ('ccp.molecule.LabeledMolecule.UniformAtomLabel',
    #                'ccp.molecule.LabeledMolecule.SingleAtomLabel',
    #                'ccp.molecule.LabeledMolecule.ResLabelFraction'):
    #    notifyFunc(self.updateSampleComponentsAfter,clazz, func)

    # Sample conditions

    for func in ('','__init__','delete','setName','setDetails','addSampleCondition','removeSampleCondition'):
      notifyFunc(self.updateSampleConditionSetsAfter,'ccp.nmr.Nmr.SampleConditionSet', func)

    for func in ('','__init__','delete','setCondition','setValue','setError','setUnit'):
      notifyFunc(self.updateSampleConditionsAfter,'ccp.nmr.Nmr.SampleCondition', func)

    # Experiments

    for func in ('','__init__','delete','setName', 'setNumDim','setSample',
                 'setSampleConditionSet', 'setSpectrometer', 'setProbe',
                 'setSampleConditionSet', 'setShiftList'):
      notifyFunc(self.updateExperimentsAfter,'ccp.nmr.Nmr.Experiment', func)

    # NMR spectrometers

    for func in ('','__init__','delete','setName',
                 'setSerialNumber','setDetails','setManufacturer',
                 'setModel','setNominalFreq',
                 'setProtonFreq','setExperiments',
                 'addExperiment','removeExperiment'):
      notifyFunc(self.updateInstrumentsAfter,'ccp.general.Instrument.NmrSpectrometer',func)

    # NMR probes

    for func in ('','__init__','delete','setName',
                 'setSerialNumber','setDetails','setManufacturer',
                 'setModel','setType', 'setDiameter', 'setExperiments',
                 'addExperiment','removeExperiment'):
      notifyFunc(self.updateInstrumentsAfter,'ccp.general.Instrument.NmrProbe',func)

    # Structures etc.

    for func in ('','__init__','delete','addEntry','removeEntry','setStructureEnsemble','setNmrConstraintStore'):
      notifyFunc(self.updateStrucGenAfter,'ccp.nmr.Nmr.StructureGeneration', func)

    for func in ('','__init__','delete','setDetails'):
      notifyFunc(self.updateStrucEnsAfter,'ccp.molecule.MolStructure.StructureEnsemble', func)

    for func in ('','__init__','delete','setName','setDetails'):
      for clazz in ('ccp.nmr.NmrConstraint.ChemShiftConstraintList',
                    'ccp.nmr.NmrConstraint.CsaConstraintList',
                    'ccp.nmr.NmrConstraint.DihedralConstraintList',
                    'ccp.nmr.NmrConstraint.DistanceConstraintList',
                    'ccp.nmr.NmrConstraint.HBondConstraintList',
                    'ccp.nmr.NmrConstraint.JCouplingConstraintList',
                    'ccp.nmr.NmrConstraint.RdcConstraintList'):
        notifyFunc(self.updateConstraintListsAfter,clazz, func)


  def selectTab(self, index):

    funcsDict = {0:(self.updateMain,),
                 1:(self.updatePeople, self.updateGroups, self.updateAddresses),
                 2:(self.updateCitations, self.updateSoftware),
                 3:(self.updateNmrLists,),
                 4:(self.updateShiftRefs, self.updateCheckShiftListPulldown),
                 5:(self.updateMolecules, self.updateChains),
                 6:(self.updateSources,),
                 7:(self.updateSamples, self.updateSampleComponents),
                 8:(self.updateSampleConditionSets, self.updateSampleConditions),
                 9:(self.updateInstruments,),
                10:(self.updateExperiments,),
                11:(self.updateStrucGens,self.updateStrucEns,self.updateConstraintSetPulldown)}

    for func in funcsDict[index]:
      func()

  def updateAll(self, project=None):
    # Called by Extend-NMR GUI to change/set CCPN projects 

    if project:
      self.project = project
      self.nmrProject = project.currentNmrProject

    self.updateEntries()
    self.checkEntryConsistency()

    self.selectTab(self.tabbedFrame.selected)    

    self.waiting = False

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    Frame.destroy(self)


  # Entry and Main tab functions

  def deleteCurrentEntry(self):

    entry = self.entry
    if entry:
      msg = 'Really delete deposition "%s" ?' % entry.name
      if showOkCancel('Confirm', msg, parent=self):
        entry.delete()
        self.changeEntry(None)

  def addNewEntry(self):

    if not self.project:
      showWarning('Failure','You need a CCPN project.', parent=self)
      return

    msg = 'Working name for deposition entry:'
    name = askString('Query', msg, parent=self) or ''
    name.strip()

    if len(name.split()) > 1:
      showWarning('Failure','Name cannot contain whitespace.', parent=self)
      return

    if name:
      store = getEntryStore(self.project)

      if store.findFirstEntry(name=name):
        showWarning('Failure','Name already in use.', parent=self)
        return

      entry = store.newEntry(name=name)
      self.changeEntry(entry)

  def changeEntry(self, entry):

    self.entry = entry
    if entry:
      self.molSystem = entry.molSystem

    else:
      self.molSystem = None

    self.updateAll()

  def checkEntryConsistency(self):

    entry = self.entry

    if entry:
      for peakList in entry.peakLists:
        experiment = peakList.dataSource.experiment

        if experiment not in entry.experiments:
          entry.addExperiment(experiment)

      # Note: can add this code to check that linked shiftLists and experiments are in the Entry.

      #for experiment in entry.experiments:
      #  shiftList = experiment.shiftList

      #  if shiftList and shiftList not in entry.findAllMeasurementLists(className="ShiftList"):
      #    entry.addMeasurementList(shiftList)

      constraintSets = set()
      for structureAnalysis in entry.structureAnalyses:
        rSet = structureAnalysis.nmrConstraintStore
        if rSet:
          constraintSets.add(rSet)

      for structureAnalysis in self.entry.structureGenerations:
        rSet = structureAnalysis.nmrConstraintStore
        if rSet:
          constraintSets.add(rSet)

      for constraintSet in constraintSets:
        for constraintList in constraintSet.constraintLists:
          for experiment in constraintList.experiments:
            if experiment not in entry.experiments:
              entry.addExperiment(experiment)

      methodStore = getMethodStore(self.project)
      entry.software = methodStore.sortedSoftware()

  def updateEntryAfter(self, entry=None):

    if self.waiting:
      return 

    else:
      self.waiting = True
      self.after_idle(self.updateAll)  

  def updateEntries(self):

    names = []
    index = []
    store = self.project.currentNmrEntryStore or \
              self.project.newNmrEntryStore(name='eciDefault')
    entries = store.sortedEntries()
    entry = self.entry

    if entries:
      names = [e.name for e in entries]
      if entry not in entries:
        entry = entries[0]
      index = entries.index(entry)  

    else:
      entry = None

    if entry is not self.entry:
      self.changeEntry(entry)

    self.entryPulldown.setup(names, entries, index)

  def updateReport(self):

    from ccpnmr.eci.CompletenessCheck import checkNmrEntryCompleteness

    textMatrix = []
    objectList = []
    colorMatrix = []

    if self.entry:
      results = checkNmrEntryCompleteness(self.entry, submissionType=self.submissionType)

      for ccpn, notes, color in results:

        datum = [ccpn, notes]
        textMatrix.append(datum)

        objectList.append(ccpn)

        colors = [None, color]
        colorMatrix.append(colors)

    self.reportTable.update(textMatrix=textMatrix,
                            objectList=objectList,
                            colorMatrix=colorMatrix)

  def updateMain(self):

    if self.entry:
      title = self.entry.title or ''
      details = self.entry.details or ''
      process = self.entry.bmrbProcessing or ''
      keywords = self.entry.keywords or []
      self.nmrMethodPulldown.set(self.entry.entryType)
    else:
      title = details = process = keywords = ''

    self.titleEntry.set(title)
    self.detailsText.setText(details)
    self.processText.setText(process)
    self.pdbKeywordsMulti.set(keywords)
    self.updateMolSystemPulldown()
    self.updateReport()

  def updateEntryTitle(self, event=None):

    text = self.titleEntry.get() or ''
    text.strip()

    if self.entry:
      self.entry.setTitle(text or None)

  def updateMolSystemPulldown(self):

    index = 0
    names = []
    molSystems = []

    if self.entry:
      molSystems = self.project.sortedMolSystems()
      molSystem = self.entry.molSystem

      if molSystems:
        names = [ms.code for ms in molSystems]

        if molSystem not in molSystems:
          molSystem = molSystems[0]

        index = molSystems.index(molSystem)  

      else:
        molSystem = None

      if molSystem is not self.entry.molSystem:
        self.changeMolSystem(molSystem)

    self.molSystemPulldown.setup(names, molSystems, index)

  def getPdbKeywords(self, entry):

    if entry:
      self.pdbKeywordsMulti.setValues(entry.keywords)

  def setEntryType(self, text):

    if self.entry:
      self.entry.entryType = text

  def setSubmissionType(self, text):

    if self.entry:
      self.submissionType = text
      self.updateReport()

  # ADIT-NMR and AutoDep have different options here -
  #   so it is a bit confusing about what to do for these 4 routines.

  def setCoordDate(self, text):

    nmrStarText = None

    if text == 'release coordinates immediatedly':
      nmrStarText = 'RELEASE NOW'
    elif text == 'release coordinates upon publication (maximum hold period 1 year)':
      nmrStarText = 'HOLD FOR PUBLICATION'
    elif text == 'place coordinates on hold for one year':
      nmrStarText = 'HOLD FOR 1 YEAR'

    if self.entry and nmrStarText:
      appData = self.entry.findFirstApplicationData(application='nmrStar',keyword='depRelCoord')
      if appData:
        self.entry.removeApplicationData(appData)

        keywds = {'application': 'nmrStar',
                  'keyword':     'depRelCoord',
                  'value':       nmrStarText}

        appData = Implementation.AppDataString(**keywds)
        self.entry.addApplicationData(appData)

  def setRestrDate(self, text):

    nmrStarText = None

    if text == 'release NMR restraints immediatedly (3 weeks after completion of deposition)':
      nmrStarText = 'RELEASE NOW'
    elif text == 'release NMR restraints data upon publication of paper':
      nmrStarText = 'HOLD FOR PUBLICATION'
    elif text == 'place one year hold on the NMR restraints':
      nmrStarText = 'HOLD FOR 1 YEAR'

    if self.entry and nmrStarText:
      appData = self.entry.findFirstApplicationData(application='nmrStar',keyword='depRelConstr')
      if appData:
        self.entry.removeApplicationData(appData)

        keywds = {'application': 'nmrStar',
                  'keyword':     'depRelConstr',
                  'value':       nmrStarText}

        appData = Implementation.AppDataString(**keywds)
        self.entry.addApplicationData(appData)

  def setNmrDate(self, text):

    nmrStarText = None

    if text == 'release other NMR data immediatedly (x weeks after completion of deposition)':
      nmrStarText = 'RELEASE NOW'
    elif text == 'release other NMR data upon publication of paper':
      nmrStarText = 'HOLD FOR PUBLICATION'
    elif text == 'place one year hold on the other NMR data':
      nmrStarText = 'HOLD FOR 1 YEAR'

    if self.entry and nmrStarText:
      appData = self.entry.findFirstApplicationData(application='nmrStar',keyword='depRelNmr')
      if appData:
        self.entry.removeApplicationData(appData)

        keywds = {'application': 'nmrStar',
                  'keyword':     'depRelNmr',
                  'value':       nmrStarText}

        appData = Implementation.AppDataString(**keywds)
        self.entry.addApplicationData(appData)

  def setSeqDate(self, text):

    nmrStarText = None

    if text == 'yes':
      nmrStarText = 'RELEASE NOW'
    elif text == 'no':
      nmrStarText = 'HOLD FOR RELEASE'

    if self.entry and nmrStarText:
      appData = self.entry.findFirstApplicationData(application='nmrStar',keyword='depRelSeq')
      if appData:
        self.entry.removeApplicationData(appData)

        keywds = {'application': 'nmrStar',
                  'keyword':     'depRelSeq',
                  'value':       nmrStarText}

        appData = Implementation.AppDataString(**keywds)
        self.entry.addApplicationData(appData)

  def setDetails(self, event=None):

    if self.entry:
      text = self.detailsText.getText()[:-1] or ''
      text.rstrip()
      if text != '':
        self.entry.details = text or None

  def setProcessing(self, event=None):

    if self.entry:
      text = self.processText.getText()[:-1] or ''
      text.rstrip()
      self.entry.bmrbProcessing = text or None

  def setPdbKeywords(self, event):

    if self.entry:
      keywords = [w.strip() for w in self.pdbKeywordsMulti.get() if w.strip()]
      self.entry.keywords = keywords

  def changeMolSystem(self, molSystem):

    if self.entry and (molSystem is not self.entry.molSystem):
      self.entry.molSystem = molSystem
      self.molSystem = molSystem

      if molSystem:

        experiments = self.entry.experiments
        for experiment in experiments:
          if molSystem not in experiment.molSystems:
            self.entry.removeExperiment(experiment)

        peakLists = self.entry.peakLists
        for peakList in peakLists:
          experiment = peakList.dataSource.experiment
          if molSystem not in experiment.molSystems:
            self.entry.removePeakList(peakList)

      self.updateChains()
      self.updateMolecules()
      self.updateCheckShiftListPulldown()

  # Should go to the right bit of the GUI according to obj type...

  def remedyReport(self, ccpnObj):

    if ccpnObj.startswith('BMRB Entry'):
      pass
    elif ccpnObj.startswith('Deposition Authors'):
      self.tabbedFrame.select(1)
    elif ccpnObj.startswith('Contact Person'):
      self.tabbedFrame.select(1)
    elif ccpnObj.startswith('Citation'):
      self.tabbedFrame.select(2)
    elif ccpnObj.startswith('Main Citation'):
      self.tabbedFrame.select(2)
    elif ccpnObj.startswith('Software'):
      self.tabbedFrame.select(2)
    elif ccpnObj.startswith('Chemical Shift Lists'):
      self.tabbedFrame.select(3)
    elif ccpnObj.startswith('Chemical Shifts'):
      self.tabbedFrame.select(3)
    elif ccpnObj.startswith('Shift List Software'):
      self.tabbedFrame.select(3)
    elif ccpnObj.startswith('Shift List Experiments'):
      self.tabbedFrame.select(10)
    elif ccpnObj.startswith('Chemical Shift References'):
      self.tabbedFrame.select(4)
    elif ccpnObj.startswith('Molecular System'):
      self.tabbedFrame.select(5)
    elif ccpnObj.startswith('Macromolecule'):
      self.tabbedFrame.select(5)
    elif ccpnObj.startswith('Biological Sources'):
      self.tabbedFrame.select(6)
    elif ccpnObj.startswith('Experimental Source'):
      self.tabbedFrame.select(6)
    elif ccpnObj.startswith('Natural Source'):
      self.tabbedFrame.select(6)
    elif ccpnObj.startswith('Sample Condition'):
      self.tabbedFrame.select(8)
    elif ccpnObj.startswith('Sample'): # Bit lazy - must be below Sample Condition
      self.tabbedFrame.select(7)
    elif ccpnObj.startswith('NMR Spectrometer'):
      self.tabbedFrame.select(9)
    elif ccpnObj.startswith('Spectrometer Manufacturer'):
      self.tabbedFrame.select(9)
    elif ccpnObj.startswith('NMR Experiment'):
      self.tabbedFrame.select(10)
    elif ccpnObj.startswith('Structural Calculation Data'):
      self.tabbedFrame.select(11)
    elif ccpnObj.startswith('Structure Generation Software'):
      self.tabbedFrame.select(11)
    elif ccpnObj.startswith('Structural Ensemble'):
      self.tabbedFrame.select(11)
    elif ccpnObj.startswith('Ensemble Models'):
      self.tabbedFrame.select(11)
    elif ccpnObj.startswith('NMR Constraint Store'):
      self.tabbedFrame.select(11)
    elif ccpnObj.startswith('Distance Constraint'):
      self.tabbedFrame.select(11)



  # NMR List functions

  def toggleNmrListSubmit(self, nmrListType):

    baseName, listName = self.classSubmitDict[self.entryObjType]

    if self.editObject in getattr(self.entry, listName):
      func = getattr(self.entry,'remove'+baseName)
    else:
      func = getattr(self.entry,'add'+baseName)

    func(self.editObject)      

  def submitAllNmrLists(self):

    baseName, listName = self.classSubmitDict[self.entryObjType]

    func = getattr(self.entry,'add'+baseName)
    for obj in getNmrLists(self.nmrProject, self.entryObjType):
      if obj not in getattr(self.entry, listName):
        func(obj)

  def setNmrListName(self, event):

    if self.editObject:

      text = self.nmrListNameEntry.get() or ''

      text.strip()

      if text:
        self.editObject.name = text

        self.updateNmrListSelectTable()

  def setNmrListSoftware(self, event):

    if self.editObject:

      software = self.nmrListSoftwarePulldown.getObject()

      if software:
        method = software.findFirstMethod()
        if not method:
          showWarning('Failure','Please set a method for this software in the References tab.', parent=self)
          self.tabbedFrame.select(2)
          return

        else:
          self.editObject.__dict__['method'] = method
          self.updateNmrListSelectTable()

      else:
        self.editObject.__dict__['method'] = None
        self.updateNmrListSelectTable()

  def getNmrListName(self, editObject):

    if editObject:
      self.nmrListNameEntry.set(editObject.name)

  def getNmrListSoftware(self, editObject):

    if editObject:

      methodStore = getMethodStore(self.project)

      index = 0
      names = []
      softs = []

      for soft in methodStore.sortedSoftware():
        names.append(soft.name)
        softs.append(soft)

      names.append('<None>')
      softs.append(None)

      if hasattr(editObject, 'method') and editObject.method:
        index = softs.index(editObject.method.software)

      self.nmrListSoftwarePulldown.setup(names, softs, index)

  def submitSelectedNmrList(self):

    baseName, listName = self.classSubmitDict[self.entryObjType]

    func = getattr(self.entry,'add'+baseName)
    if self.editObject and self.editObject not in getattr(self.entry, listName):
      func(self.editObject)

  def withdrawSelectedNmrList(self):

    baseName, listName = self.classSubmitDict[self.entryObjType]

    func = getattr(self.entry,'remove'+baseName)
    if self.editObject and self.editObject in getattr(self.entry, listName):
      func(self.editObject)

  def withdrawAllNmrLists(self):

    baseName, listName = self.classSubmitDict[self.entryObjType]

    func = getattr(self.entry,'remove'+baseName)
    for obj in getNmrLists(self.nmrProject, self.entryObjType):
      if obj in getattr(self.entry, listName):
        func(obj)

  def selectGenTable(self, obj, row, col):

    self.entryObjType = obj
    self.updateNmrListSelectTable()

  def doEditTableEditMarkExtraRules(self, obj, row, col):
      
    if (col == 4) and (obj.className == 'PeakList'):
      return False

    return True

  def selectEditTable(self, obj, row, col):

    self.editObject = obj

  def updateNmrLists(self, entry=None):  

    if entry and (entry is not self.entry):
      return

    headColor = [None] * 4

    textMatrix = []
    colorMatrix = []
    objectList = []

    if self.entry:
      for className, title, description in NMR_LISTS_DATA:
        selected = getattr(self.entry, self.classSubmitDict[className][1])
        selected = [s for s in selected if s.className == className]
        nmrLists = getNmrLists(self.nmrProject, className)

        textMatrix.append([title,len(nmrLists),len(selected),description])
        colorMatrix.append(headColor)
        objectList.append(className)

    self.nmrListsTable.update(textMatrix=textMatrix,
                              colorMatrix=colorMatrix,
                              objectList=objectList)

    self.updateNmrListSelectTable()

  def updateNmrListSelectTable(self):

    textMatrix = []
    objectList = []
    headingList = []
    colorMatrix = []
    baseName, listName = self.classSubmitDict[self.entryObjType]

    if self.entry:
      for obj in getNmrLists(self.nmrProject, self.entryObjType):
        if obj in getattr(self.entry, listName):
          use = 'Yes'
          colors = [NICE_GREEN] * 5
        else:
          use = 'No'
          colors = [None] * 5

        if hasattr(obj, 'measurements'):
          size = len(obj.measurements)
        elif hasattr(obj, 'peaks'):
          size = len(obj.peaks)
        elif hasattr(obj, 'derivations'):
          size = len(obj.derivations)
        else:
          size = None

        softName = None
        method = None

        if obj and hasattr(obj, 'method'): # Peak lists don't have a method attribute
          method = obj.method

        if method:
          if method.software:
            softName = method.software.name

        datum = [getNmrListName(obj),
                 obj.name,
                 size,
                 use,
                 softName]

        objectList.append(obj)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    tipTexts = ['Indicate the type and name of each list, based on list type selectd in the top frame',
                'Name of list which can also be changed',
                'Number of data points in the list',
                'Is the list selected for deposition? - toggle by clicking',
                'Software used to make the list - the software will also need a method selected in the References tab'
                ]

    headingList = self.nmrListSelectTable.headingList
    headingList[0] = self.entryObjType

    self.nmrListSelectTable.update(textMatrix=textMatrix,
                                   objectList=objectList,
                                   colorMatrix=colorMatrix,
                                   headingList=headingList,
                                   tipTexts=tipTexts)

  def addAllNmrLists(self):

    if self.entry:
      for className, title, description in NMR_LISTS_DATA:
         baseName, listName = self.classSubmitDict[className]

         func = getattr(self.entry,'add'+baseName)
         for obj in getNmrLists(self.nmrProject, className):
           if obj not in getattr(self.entry, listName):
             func(obj)

  def removeAllNmrLists(self):

    if self.entry:
      for className, title, description in NMR_LISTS_DATA:
        baseName, listName = self.classSubmitDict[className]

        func = getattr(self.entry,'remove'+baseName)
        for obj in getNmrLists(self.nmrProject, className):
          if obj in getattr(self.entry, listName):
            func(obj)


  # People functions

  def updatePersonAfter(self, person):

    affStore = getAffiliationStore(self.project)

    if affStore and self.entry:
      if person in affStore.persons:
        self.updatePeople()

      if person in self.entry.contactPersons:
        self.updateGroups()
        self.updateAddresses()

  def updatePeople(self):

    #if self.entry:
    #  if self.author not in self.entry.authors:
    #    self.author = None 
    #else:
    #  self.author = None

    textMatrix = []
    objectList = []

    affStore = getAffiliationStore(self.project)

    if self.entry and affStore:
      #for i, author in enumerate(self.entry.authors):
      for i, person in enumerate(affStore.sortedPersons() ):
        gpName = None

        cPinG = getPersonInGroup(person)

        if cPinG:
          gpName = cPinG.group.name

        appData = person.findFirstApplicationData(application='AutoDep',keyword='Primary')

        datum = [i+1,
                 person in self.entry.authors and 'Yes' or None,
                 person in self.entry.contactPersons and 'Yes' or None,
                 appData and 'Yes' or None,
                 person.familyName,
                 person.title,
                 person.givenName,
                 ' '.join(person.middleInitials),
                 person.familyTitle,
                 gpName]

        objectList.append(person)
        textMatrix.append(datum)

    self.personTable.update(textMatrix=textMatrix,
                            objectList=objectList)

  def selectPersonTable(self, obj, row, col):

    self.person = obj

  def addPerson(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    fName = askString('Query','Family name:',parent=self) or ''
    fName.strip()

    if fName:
      gName = askString('Query','Given name:',parent=self) or ''
      gName.strip()

      if gName:
        affStore = getAffiliationStore(self.project)
        person = affStore.findFirstPerson(familyName=fName, givenName=gName)
        if not person: # or (person in self.entry.authors):
          person = affStore.newPerson(familyName=fName, givenName=gName) #, title=TITLES[0])

      #self.entry.addAuthor(person)

  def removePerson(self):

    affStore = getAffiliationStore(self.project)

    if not self.entry and affStore:
      return

    person = self.person
    if person:

      msg = 'Delete person entirely for %s?'
      if showYesNo('Query', msg % getPersonName(person), parent=self):
        person.delete()
        self.updatePeople()

  def setPersonAsContact(self):

    if not (self.entry and self.person):
      return

    authors = list(self.entry.authors)

    if self.person not in authors:
      showWarning('Failure','Please set person as an author.', parent=self)
      return

    contactPersons = self.entry.sortedContactPersons()

    if not self.person in contactPersons:
      contactPersons.append(self.person)

      self.entry.contactPersons = contactPersons

  def setPersonAsAuthor(self):

    if not (self.entry and self.person):
      return

    authors = list(self.entry.authors)

    if not self.person in authors:
      authors.append(self.person)

      self.entry.authors = authors

  def setPersonAsPrimary(self):

    if not (self.entry and self.person):
      return

    contactPersons = self.entry.sortedContactPersons()

    if self.person not in contactPersons:
      showWarning('Failure','Please set person as a contact.', parent=self)
      return

    for cp in contactPersons:
      appData = cp.findFirstApplicationData(application='AutoDep',keyword='Primary')
      if appData:
        cp.removeApplicationData(appData)

    keywds = {'application': 'AutoDep',
              'keyword':     'Primary',
              'value':       'Yes'}

    appData = Implementation.AppDataString(**keywds)
    self.person.addApplicationData(appData)

  def setFamilyName(self, event):

    if not (self.entry and self.person):
      return

    text = self.familyNameEntry.get() or ''
    text.strip()

    if text:
      self.person.familyName = text

  def setPersonTitle(self, event):

    if not (self.entry and self.person):
      return

    title = self.titlePulldown.getObject()
    self.person.title = title

  def setGivenName(self, event):

    if not (self.entry and self.person):
      return

    text = self.givenNameEntry.get() or ''
    text.strip()

    if text:
      self.person.givenName = text

  def setInitials(self, event):

    if not (self.entry and self.person):
      return

    text = self.initialsEntry.get() or ''
    initials = []

    for letter in text:
      if letter.upper() != letter.lower():
        initials.append(letter.upper())

    self.person.setMiddleInitials(initials)

  def setFamilyTitle(self, event):

    if not (self.entry and self.person):
      return

    title = self.famTitlePulldown.getObject()
    self.person.familyTitle = title

  def setPersonGroup(self, event):

    if not (self.entry and self.person):
      return

    group = self.personGroupPulldown.getObject()

    personInGroup = setPersonInGroup(self.person, group=group)

  def getFamilyName(self, person):

    self.familyNameEntry.set(person.familyName)

  def getPersonTitle(self, person):

    title = person.title

    if not title:
      return

    names = TITLES + ['None>',]
    objects = TITLES + [None,]

    if title not in objects:
      names.append(title)
      objects.append(title)

    index = names.index(title)
    self.titlePulldown.setup(names, objects, index)  

  def getGivenName(self, person):

    self.givenNameEntry.set(person.givenName)

  def getInitials(self, person):

    text = ' '.join(person.middleInitials)
    self.initialsEntry.set(text)

  def getFamilyTitle(self, person):

    title = person.familyTitle

    if not title:
      return

    names = FAM_TITLES + ['<None>',]
    objects = FAM_TITLES + [None,]

    if title not in objects:
      names.append(title)
      objects.append(title)

    index = objects.index(title)
    self.famTitlePulldown.setup(names, objects, index)  

  def getPersonGroup(self, person):

    if not self.entry:
      return

    names = []
    groups = []
    index = 0

    for group in self.entry.sortedLaboratories():
      names.append(group.name)
      groups.append(group)

    if self.group:
      if self.group not in groups:
        groups.insert(0, self.group)
        names.insert(0, self.group.name)

      index = groups.index(self.group)

    groups.append(None)
    names.append('<None>')

    self.personGroupPulldown.setup(names, groups, index)

  def toggleContact(self, person):

    authors = list(self.entry.authors)

    if person not in authors:
      showWarning('Failure','Please set person as an author.', parent=self)
      return

    contactPersons = self.entry.sortedContactPersons()

    if person in contactPersons:
      contactPersons.remove(person)
      appData = person.findFirstApplicationData(application='AutoDep',keyword='Primary')
      if appData:
        person.removeApplicationData(appData)

    else:
      contactPersons.append(person)

    self.entry.contactPersons = contactPersons

  def toggleAuthor(self, person):

    authors = list(self.entry.authors)
    contactPersons = self.entry.sortedContactPersons()

    if person in authors:
      authors.remove(person)
      if person in contactPersons:
        contactPersons.remove(person)
        self.entry.contactPersons = contactPersons
      appData = person.findFirstApplicationData(application='AutoDep',keyword='Primary')
      if appData:
        person.removeApplicationData(appData)

    else:
      authors.append(person)

    self.entry.authors = authors

  def togglePrimary(self, person):

    contactPersons = self.entry.sortedContactPersons()

    if person not in contactPersons:
      showWarning('Failure','Please set person as a contact.', parent=self)
      return

    primaryFlag = False

    appData = person.findFirstApplicationData(application='AutoDep',keyword='Primary')

    if appData:
      primaryFlag = True

    for cp in contactPersons:
      appData = cp.findFirstApplicationData(application='AutoDep',keyword='Primary')
      if appData:
        cp.removeApplicationData(appData)

    if not primaryFlag:
      keywds = {'application': 'AutoDep',
                'keyword':     'Primary',
                'value':       'Yes'}

      appData = Implementation.AppDataString(**keywds)
      person.addApplicationData(appData)


  # Group functions

  def updateGroupAfter(self, group):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateGroups)

  def updateGroups(self):

    entry = self.entry

    if entry:
      if self.group not in self.entry.laboratories:
        self.group = None
    else:
      self.group = None

    textMatrix = []
    objectList = []

    if entry:
      for i, group in enumerate(entry.sortedLaboratories() ):
        organisation = group.organisation

        stateProv = None

        appData = organisation.findFirstApplicationData(application='nmrStar',keyword='State')

        if appData:
          stateProv = appData.value

        addresses = ', '.join(organisation.addresses)
        city = organisation.city
        country = organisation.country
        postalCode = organisation.postalCode
        orgType = organisation.organisationType

        datum = [group.name,
                 addresses,
                 city,
                 stateProv,
                 postalCode,
                 country,
                 orgType]

        objectList.append(group)
        textMatrix.append(datum)

    self.groupTable.update(textMatrix=textMatrix,
                           objectList=objectList)

    self.waiting = False

  def selectGroupTable(self, obj, row, col):

    self.group = obj

  def addGroup(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    affStore = getAffiliationStore(self.project)

    gpName = askString('Query','Group name:',parent=self) or ''
    gpName.strip()

    if gpName:
      organisation = affStore.findFirstOrganisation(name=gpName) or \
                     affStore.newOrganisation(name=gpName)

      group = organisation.findFirstGroup(name=gpName) or \
              organisation.newGroup(name=gpName)

      self.entry.addLaboratory(group)

  def removeGroup(self):

    if not self.entry:
      return

    group = self.group
    
    if group:
      self.entry.removeLaboratory(group)
      msg = 'Delete group for %s entirely?'
      if showYesNo('Query', msg % group.name, parent=self):
        group.delete()
        group = None

      affStore = getAffiliationStore(self.project)
      for person in affStore.sortedPersons():
        pinG = getPersonInGroup(person)
        if pinG:
          person.setCurrentPersonInGroup(None)
          pinG.delete()
          ping = None

  def getMailAddress(self, group):

    if group:
      self.mailAddressEntry.set(','.join(group.organisation.addresses) )

  def getGroupCity(self, group):

    if group:
      self.groupCityEntry.set(group.organisation.city)

  def getStateProv(self, group):
  
    if group:
      orgn = group.organisation

      value = None

      appData = orgn.findFirstApplicationData(application='nmrStar',keyword='State')

      if appData:
        value = appData.value

      if not value:
        return

      self.stateProvEntry.set(value)

  def getCountry(self, group):

    if group:
      self.countryEntry.set(group.organisation.country)

  def getPostalCode(self, group):

    if group:
      self.postalCodeEntry.set(group.organisation.postalCode)

  def getOrgType(self, group):

    orgType = group.organisation.organisationType

    if not orgType:
      return

    names = ORGN_TYPES + ['<None>',]
    objects = ORGN_TYPES + [None,]

    if orgType not in objects:
      names.append(orgType)
      objects.append(orgType)

    index = objects.index(orgType)
    self.orgTypePulldown.setup(names, objects, index)

  def setMailAddress(self, event):

    if self.group:
      text = self.mailAddressEntry.get() or ''
      text.strip()

      if text:
        self.group.organisation.addresses = text.split(',')

  def setGroupCity(self, event):

    if self.group:
      text = self.groupCityEntry.get() or ''
      text.strip()

      if text:
        self.group.organisation.city = text

  def setStateProv(self, event):

    if self.group:
      orgn = self.group.organisation

      appData = orgn.findFirstApplicationData(application='nmrStar',keyword='State')

      if appData:
        for appData in orgn.findAllApplicationData(application='nmrStar',keyword='State'):
          orgn.removeApplicationData(appData)

      text = self.stateProvEntry.get() or ''
      text.strip()

      if text:
        keywds = {'application': 'nmrStar',
                  'keyword':     'State',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        orgn.addApplicationData(appData)

  def setCountry(self, event):

    if self.group:
      text = self.countryEntry.get() or ''
      text.strip()

      if text:
        self.group.organisation.country = text

  def setPostalCode(self, event):

    if self.group:
      text = self.postalCodeEntry.get() or ''
      text.strip()

      if text:
        self.group.organisation.postalCode = text

  def setOrgType(self, event):

    orgType = self.orgTypePulldown.getObject()
    self.group.organisation.organisationType = orgType


  # Contact person address functions

  def updateAddressAfter(self, person):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateAddresses)

  def updateAddresses(self):

    entry = self.entry

    if entry:
      if self.contactPerson not in self.entry.contactPersons:
        self.contactPerson = None
    else:
      self.contactPerson = None

    textMatrix = []
    objectList = []

    if entry:
      for i, contactPerson in enumerate(entry.sortedContactPersons() ):
        email = phoneNumber = faxNumber = position = None

        cPinG = getPersonInGroup(contactPerson)
        
        if cPinG:
          email = cPinG.emailAddress
          phoneNumber = '/'.join(cPinG.phoneNumbers)
          faxNumber = cPinG.faxNumber
          position = cPinG.position

        datum = [contactPerson.familyName,
                 contactPerson.givenName,
                 email,
                 phoneNumber,
                 faxNumber,
                 position]

        objectList.append(contactPerson)
        textMatrix.append(datum)

    self.addressTable.update(textMatrix=textMatrix,
                             objectList=objectList)

    self.waiting = False

  def selectAddressTable(self, obj, row, col):

    self.contactPerson = obj

  def addAddress(self):

    if not self.entry or not self.contactPerson:
      showWarning('Failure','Please select a person from the list.', parent=self)
      return

    person = self.contactPerson
    personInGroup = getPersonInGroup(person)

    if not personInGroup:
      showWarning('Failure','Please ensure the person belongs to a group in the other frames.', parent=self)
      return

    while True:
      email = askString('Query','Email address:',parent=self) or ''
      email.strip()
      if email.count('@') != 1: # or not email.count('.'): # Do email addresses need at least one '.'
        showWarning('Failure','The email address is of the form "abc@this.address".', parent=self)
        continue
      break

    if email and personInGroup:
      personInGroup.emailAddress = email

  def removeAddress(self):

    if not self.entry or not self.contactPerson:
      return

    person = self.contactPerson
    personInGroup = person.currentPersonInGroup
    
    if person and personInGroup:
      msg = 'Delete address for %s entirely?'
      if showYesNo('Query', msg % getPersonName(person), parent=self):
        person.setCurrentPersonInGroup(None)
        personInGroup.delete()
        self.contactPerson = None

  def doAddressTableEditMarkExtraRules(self, obj, row, col):

    if not self.entry or not self.contactPerson:
      return False

    personInGroup = getPersonInGroup(self.contactPerson)

    if not personInGroup:
      return False

    return True

  def setEmailAddress(self, event):

    person = self.contactPerson
    if person:
      text = self.emailAddressEntry.get() or ''
      text.strip()

      if text:
        personInGroup = getPersonInGroup(person)
        personInGroup.emailAddress = text

  def setTelephone(self, event):

    if self.contactPerson:
      personInGroup = getPersonInGroup(self.contactPerson)
      if personInGroup is None or not personInGroup.emailAddress:
        showWarning('Failure','Please add an email address using the button below first.', parent=self)
        return
        
      text = self.telephoneEntry.get() or ''
      text.strip()

      if text:
        personInGroup.phoneNumbers = [text,]

  def setFax(self, event):

    if self.contactPerson:
      personInGroup = getPersonInGroup(self.contactPerson)
      if personInGroup is None or not personInGroup.emailAddress:
        showWarning('Failure','Please add an email address using the button below first.', parent=self)
        return
        
      text = self.faxEntry.get() or ''
      text.strip()

      if text:
        personInGroup.faxNumber = text

  def setPosition(self, event):

    if self.contactPerson:
      personInGroup = getPersonInGroup(self.contactPerson)
      if personInGroup is None or not personInGroup.emailAddress:
        showWarning('Failure','Please add an email address using the button below first.', parent=self)
        return
        
      orgType = self.positionPulldown.getObject()
      personInGroup = getPersonInGroup(self.contactPerson) 
      personInGroup.position = orgType

  def getEmailAddress(self, person):

    if person and person.currentPersonInGroup:
      self.emailAddressEntry.set(person.currentPersonInGroup.emailAddress)

  def getTelephone(self, person):

    if person and person.currentPersonInGroup:
      numbers = person.currentPersonInGroup.phoneNumbers
    
      if numbers:
        number = numbers[0]
      else:
        number = None  
        
      self.telephoneEntry.set(number)

  def getFax(self, person):

    if person and person.currentPersonInGroup:
      self.faxEntry.set(person.currentPersonInGroup.faxNumber)

  def getPosition(self, person):

    if person and person.currentPersonInGroup:

      posRole = person.currentPersonInGroup.position

      if not posRole:
        return

      names = POS_ROLES + ['<None>',]
      objects = POS_ROLES + [None,]

      if posRole not in objects:
        names.append(posRole)
        objects.append(posRole)

      index = objects.index(posRole)
      self.positionPulldown.setup(names, objects, index)


  # Chem shift reference functions

  def updateShiftRefs(self, obj=None):

    objectList  = []
    textMatrix  = []
    if self.nmrProject:
      for shiftReference in self.nmrProject.sortedShiftReferences():

        refClass = shiftReference.className[:8]

        if refClass == 'External':
          geometry = shiftReference.sampleGeometry
          location = shiftReference.location
          axis     = shiftReference.axis
        else:
          geometry = location = axis = None

        #' '.join([e.name for e in shiftReference.experiments]),
        datum = [shiftReference.serial,
                 refClass,
                 shiftReference.isotopeCode,
                 len(shiftReference.experiments),
                 shiftReference.molName,
                 shiftReference.atomGroup,
                 shiftReference.value,
                 shiftReference.unit,
                 shiftReference.referenceType,
                 shiftReference.indirectShiftRatio,
                 geometry,location,axis]

        textMatrix.append(datum)
        objectList.append(shiftReference)

    self.shiftRefTable.update(objectList=objectList, textMatrix=textMatrix)

    self.waiting = False

  def doShiftRefEditMarkExtraRules(self, obj, row, col):

    if (col > 9) and (obj.className != 'ExternalShiftReference'):
      return False

    return True  

  def toggleShiftRefType(self, shiftReference):

    if shiftReference:
      if shiftReference.referenceType == 'direct':
        shiftReference.setReferenceType('indirect')
      else:
        shiftReference.setReferenceType('direct')

  def getShiftRefValue(self, shiftReference):

    value = 0.0
    if shiftReference:
      value = shiftReference.value

    self.valueEntry.set(value)

  def getShiftRefRatio(self, shiftReference):

    value = 0.0
    if shiftReference:
      value = shiftReference.indirectShiftRatio

    self.ratioEntry.set(value)

  def getShiftRefGeometry(self, shiftReference):

    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.sampleGeometry

    self.geometryEntry.set(text)

  def getShiftRefLocation(self, shiftReference):

    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.location

    self.locationEntry.set(text)

  def getShiftRefAxis(self, shiftReference):

    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.axis

    self.axisEntry.set(text)

  def getShiftRefAtomGroup(self, shiftReference):

    text = ''
    if shiftReference:
      text = shiftReference.atomGroup

    self.atomGroupEntry.set(text)

  def getShiftRefIsotope(self, shiftReference):

    self.isotopePulldown.set(shiftReference.isotopeCode)

  def getShiftRefMolName(self, shiftReference):

    molName = shiftReference.molName

    if not molName:
      return

    names = SHIFT_REF_COMPOUNDS + ['<Other>',]
    objects = SHIFT_REF_COMPOUNDS + [None,]

    if molName not in objects:
      names.append(molName)
      objects.append(molName)

    index = names.index(molName)
    self.molNamePulldown.setup(names, objects, index)  

  def getShiftRefUnit(self, shiftReference):

    self.unitPulldown.set(shiftReference.unit)

  def setShiftRefValue(self, event):

    if self.shiftReference:
      value = self.valueEntry.get() or 0.0

      self.shiftReference.value = value

  def setShiftRefRatio(self, event):

    if self.shiftReference:
      value = self.ratioEntry.get() or None

      self.shiftReference.indirectShiftRatio = value

  def setShiftRefGeometry(self, event):

    if self.shiftReference:
      text = self.geometryEntry.get() or None

      self.shiftReference.sampleGeometry = text

  def setShiftRefLocation(self, event):

    if self.shiftReference:
      text = self.locationEntry.get() or None

      self.shiftReference.location = text

  def setShiftRefAxis(self, event):

    if self.shiftReference:
      text = self.axisEntry.get() or None

      self.shiftReference.axis = text

  def setShiftRefAtomGroup(self, event):

    if self.shiftReference:
      text = self.atomGroupEntry.get() or None

      self.shiftReference.atomGroup = text

  def setShiftRefIsotope(self, null):

    isotopeCode = self.isotopePulldown.getText()

    self.shiftReference.isotopeCode = isotopeCode

  def setShiftRefMolName(self, null):

    molName = self.molNamePulldown.getObject() 

    if molName is None:
      msg = 'Molecule name:'
      molName = askString('Input', msg, parent=self) or ''
      molName.strip()
        
    if not molName:
      showWarning('Failure','Name cannot be an empty string.', parent=self)
      return

    self.shiftReference.molName = molName

  def setShiftRefUnit(self, null):

    unit = self.unitPulldown.getText()

    self.shiftReference.unit = unit

  def addStdIupacRef(self):

    if self.nmrProject:

      findFirstRef = self.nmrProject.findFirstShiftReference

      newIntRef = self.nmrProject.newInternalShiftReference
      newExtRef = self.nmrProject.newExternalShiftReference

      keywds = {'isotopeCode':        '1H',
                'molName':            'DSS',
                'atomGroup':          'methyl protons',
                'unit':               'ppm',
                'value':              0.00,
                'referenceType':      'direct',
                'indirectShiftRatio': 1.000000000}

      shiftRef1H = findFirstRef(className='InternalShiftReference', **keywds)

      if not shiftRef1H:
        shiftRef1H = newIntRef(**keywds)

      self.shiftReference = shiftRef1H

      # Note we set 13C, 15N and 31P to be of class ExternalShiftReference
      # Then set the Application data to be 'n/a' so that they get exported correctly to NmrStar.

      keywds['isotopeCode']        = '13C'
      keywds['referenceType']      = 'indirect'
      keywds['indirectShiftRatio'] = 0.251449530

      appKeywds = {'application': 'nmrStar',
                   'keyword':     'refType',
                   'value':       'n/a'}

      shiftRef13C = findFirstRef(className='ExternalShiftReference', **keywds)

      if not shiftRef13C or not shiftRef13C.findFirstApplicationData(**appKeywds):
        shiftRef13C = newExtRef(**keywds)
        shiftRef13C.addApplicationData(Implementation.AppDataString(**appKeywds) )

      keywds['isotopeCode']        = '15N'
      keywds['indirectShiftRatio'] = 0.101329118

      shiftRef15N = findFirstRef(className='ExternalShiftReference', **keywds)

      if not shiftRef15N or not shiftRef15N.findFirstApplicationData(**appKeywds):
        shiftRef15N = newExtRef(**keywds)
        shiftRef15N.addApplicationData(Implementation.AppDataString(**appKeywds) )

      keywds['isotopeCode']        = '31P'
      keywds['indirectShiftRatio'] = 0.404808636

      shiftRef31P = findFirstRef(className='ExternalShiftReference', **keywds)

      if not shiftRef31P or not shiftRef31P.findFirstApplicationData(**appKeywds):
        shiftRef31P = newExtRef(**keywds)
        shiftRef31P.addApplicationData(Implementation.AppDataString(**appKeywds) )

  def addInternalShiftRef(self):

    if self.nmrProject:
      newRef = self.nmrProject.newInternalShiftReference
      self.shiftReference = newRef(isotopeCode='1H', molName='TSP',
                                   atomGroup='H', value=0.000,
                                   referenceType='direct')

  def addExternalShiftRef(self):

    if self.nmrProject:
      newRef = self.nmrProject.newExternalShiftReference
      self.shiftReference = newRef(isotopeCode='1H', molName='TSP',
                                   atomGroup='H', value=0.000,
                                   referenceType='direct')

  def removeShiftRefs(self):

    haveExpts = False
    for shiftReference in self.shiftRefTable.currentObjects:
      if shiftReference.experiments:
        haveExpts = True
        break

    msg = 'Really delete shift references with links to experiments?'
    if haveExpts and not showOkCancel('Confirm', msg, parent=self):
      return

    for shiftReference in self.shiftRefTable.currentObjects:
      shiftReference.delete()

  def selectShiftRef(self, shiftRef, row, col):

    if shiftRef:
      self.shiftReference = shiftRef


  # Shift Check functions

  def updateCheckShiftListPulldown(self):

    names = []
    index = 0

    if self.nmrProject:
      shiftLists = getNmrLists(self.nmrProject, 'ShiftList')
    else:
      shiftLists = []  

    shiftList = self.checkShiftList

    if shiftLists:
      names = [ '%s : %d' % (sl.name, sl.serial) for sl in shiftLists ]

      if shiftList not in shiftLists:
        shiftList = shiftLists[0]

      index = shiftLists.index(shiftList)  

    else:
      shiftList = None

    if shiftList is not self.checkShiftList:
      self.changeCheckShiftList(shiftList)

    self.checkShiftListPulldown.setup(names, shiftLists, index)

  def changeCheckShiftList(self, shiftList):

    if shiftList is not self.checkShiftList:
      self.checkShiftList = shiftList
      self.updateShiftCheck()

  def updateShiftCheck(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if self.checkShiftList:

      data = analyseChemicalShifts(self.checkShiftList)

      for shift, datum in data:
        #num, isotope, name, boundResonances, mean, \
        # coil, prob, value, sd, peakDelta, \
        # nContribs, isDuplicate, boundWarn = datum

        num, isotope, name, mean, coil, dev, value, sd = datum

        resonanceSet = shift.resonance.resonanceSet
        if resonanceSet:
          atom = resonanceSet.findFirstAtomSet().findFirstAtom()
          if atom.topObject is not self.molSystem:
            continue

        color = None             
        if dev is not None:
          if dev >= 3.0:
            color = '#FF0000'
          elif dev >= 2.0:
            color = '#D0D060'
          #elif dev < 0.05:
          #  color = '#d0d0a0'

        if color is None:
          continue

        colors = [color] * 7

        objectList.append(shift)
        textMatrix.append([num, isotope, name, value, sd, mean, coil])
        colorMatrix.append(colors)

    self.shiftCheckTable.update(textMatrix=textMatrix,
                                objectList=objectList,
                                colorMatrix=colorMatrix)


  # Molecule functions

  def selectChain(self, chain, row, col):

    self.chain = chain

  def getChainConfIsomer(self, chain):

    if len(chain.residues) > 1:
      names = BOOL[:]
      confIsomers = BOOL[:]
    
    else:
      names = []
      confIsomers = []
   
    names.append('undefined')
    confIsomers.append(None)

    confIsomer = chain.conformationalIsomer
    if confIsomer not in confIsomers:
       names.insert(0, confIsomer)
       confIsomers.insert(0, confIsomer)

    index = confIsomers.index(confIsomer)

    self.chainConfIsomerPulldown.setup(names, confIsomers, index)

  def setChainConfIsomer(self, event):

    if self.chain:
      confIsomer = self.chainConfIsomerPulldown.getObject() or None
      self.chain.conformationalIsomer = confIsomer

  def getChainPhysicalState(self, chain):
    
    if len(chain.residues) > 1:
      names  = POLYMER_PHYSICAL_STATES[:]
      states = POLYMER_PHYSICAL_STATES[:]
    
    else:
      names = []
      states = []
   
    names.append('undefined')
    states.append(None)
    
    state = chain.physicalState
    if state not in states:
       names.insert(0, state)
       states.insert(0, state)

    index = states.index(state)

    self.chainPhysicalStatePulldown.setup(names, states, index)

  def setChainPhysicalState(self, event):

    if self.chain:
      state = self.chainPhysicalStatePulldown.getObject() or None
      self.chain.physicalState = state

  def getChainChemExchState(self, chain):

    if len(chain.residues) > 1:
      names = BOOL[:]
      states = BOOL[:]
    
    else:
      names = []
      states = []

    names.append('undefined')
    states.append(None)
    
    state = chain.chemExchangeState
    if state not in states:
       names.insert(0, state)
       states.insert(0, state)

    index = states.index(state)

    self.chainChemExchStatePulldown.setup(names, states, index)

  def setChainChemExchState(self, event):

    if self.chain:
      state = self.chainChemExchStatePulldown.getObject() or None
      self.chain.chemExchangeState = state

  def selectMolecule(self, obj, row, col):

    self.molecule = obj

  def getMoleculeDbName(self, molecule):

    if not molecule:
      return

    if len(molecule.molResidues) > 1:
      names = DB_NAMES[:]
      dbs   = DB_NAMES[:]
    
    else:
      names = []
      dbs   = []
   
    names.append('<Other>')
    dbs.append(None)

    dbName = None

    self.align = molecule.findFirstAlignment()

    if self.align:
      dbRef = self.align.dbRef
      if dbRef:
        self.dbName = dbRef.database.name

    if self.dbName is not None and self.dbName not in dbs:
       dbs.insert(0, self.dbName)
       names.insert(0, self.dbName)

    index = dbs.index(self.dbName)

    self.moleculeDbNamePulldown.setup(names, dbs, index)

  def getMoleculeDbRef(self, molecule):

    if not molecule:
      return

    self.align = molecule.findFirstAlignment()

    if self.align:
      dbRef = self.align.dbRef

      if dbRef:
        if not dbRef.code:
          self.accNum = dbRef.name
        else:
          self.accNum = dbRef.code

        if self.accNum:
          self.moleculeDbRefEntry.set(self.accNum)

  def getMoleculeDbRefDetails(self, molecule):

    if not molecule:
      return

    self.align = molecule.findFirstAlignment()

    if self.align:
      dbRef = self.align.dbRef

      if dbRef:
        details = dbRef.details

        if details:
          self.moleculeDbRefDetailsEntry.set(details)

  def getMoleculeDomain(self, molecule):

    if not molecule:
      return

    appData = molecule.findFirstApplicationData(application='nmrStar',keyword='Fragment')

    if appData:
      value = appData.value

      self.moleculeDomainEntry.set(value)

  def getMoleculeECnum(self, molecule):

    if not molecule:
      return

    ecnums = []

    msn = molecule.findFirstMoleculeSysName(namingSystem='EC')

    if msn:
      ecnums = msn.name.split(',')

    self.moleculeECNumMulti.setValues(ecnums)

  def getMoleculeMutation(self, molecule):

    if not molecule:
      return

    appData = molecule.findFirstApplicationData(application='nmrStar',keyword='Mutation')

    if appData:
      value = appData.value

      self.moleculeMutEntry.set(value)

  def getMoleculeDetails(self, molecule):

    if not molecule:
      return

    self.moleculeDetailsEntry.set(molecule.details)

  def setMoleculeDomain(self, event):

    if self.molecule:
      text = self.moleculeDomainEntry.get() or ''
      text.strip()

      appData = self.molecule.findFirstApplicationData(application='nmrStar',keyword='Fragment')

      if appData:
        for appData in self.molecule.findAllApplicationData(application='nmrStar',keyword='Fragment'):
          self.molecule.removeApplicationData(appData)

      if text != '':

        keywds = {'application': 'nmrStar',
                  'keyword':     'Fragment',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.molecule.addApplicationData(appData)

  def setMoleculeECnum(self, event):

    if self.molecule:
      ecnumList = [w.strip() for w in self.moleculeECNumMulti.get() if w.strip()]

      ecnums = ','.join(ecnumList)

      msn = self.molecule.findFirstMoleculeSysName(namingSystem='EC')

      if not msn:
        if ecnums != '':
          self.molecule.newMoleculeSysName(name=ecnums, namingSystem='EC')

      else:
        del(self.molecule.__dict__['moleculeSysNames']['EC'])

        if ecnums != '':
          self.molecule.newMoleculeSysName(name=ecnums, namingSystem='EC')

  def setMoleculeMutation(self, event):

    if self.molecule:
      text = self.moleculeMutEntry.get() or ''
      text.strip()

      appData = self.molecule.findFirstApplicationData(application='nmrStar',keyword='Mutation')

      if appData:
        for appData in self.molecule.findAllApplicationData(application='nmrStar',keyword='Mutation'):
          self.molecule.removeApplicationData(appData)

      if text != '':

        keywds = {'application': 'nmrStar',
                  'keyword':     'Mutation',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.molecule.addApplicationData(appData)

  def setMoleculeDetails(self, event):

    if self.molecule:
      details = self.moleculeDetailsEntry.get() or ''
      details.strip()

      if details != '':
        self.molecule.details = details

  def setMoleculeDbName(self, event):

    if self.molecule:
      self.dbName = self.moleculeDbNamePulldown.getObject() or None

      if self.dbName is None:
        msg = 'Database name:'
        self.dbName = askString('Input', msg, parent=self) or ''
        self.dbName.strip()
        
      if not self.dbName:
        showWarning('Failure','Name cannot be an empty string.', parent=self)
        return

      db = self.project.findFirstDatabase(name=self.dbName)
      if not db:
        db = self.project.newDatabase(name=self.dbName)

      self.accNum = 'TmpAcc'

      #print 'DB: [%s]' % db.name

      dbRef = db.findFirstEntry(name=self.accNum)

      if not dbRef:
        dbRef = db.newEntry(name=self.accNum, code=self.accNum)
      else:
        for i in range(1,100):
          accNum = self.accNum + str(i)
          dbRef = db.findFirstEntry(name=accNum)
          if dbRef:
            continue
          
          dbRef = db.newEntry(name=accNum, code=self.accNum)
          break

      if not dbRef.code:
        dbRef.code = self.accNum

      #print 'DBREF: [%s]' % dbRef.name

      self.accNum = dbRef.name

      self.align = self.molecule.findFirstAlignment()

      if not self.align:
        self.align = self.molecule.newAlignment(dbRef=dbRef)
      else:
        self.align.__dict__['dbRef'] = dbRef
        #self.align.dbRef.delete()
        #self.align.setDbRef(dbRef)

  def setMoleculeDbRef(self, event):

    if self.molecule:
      text = self.moleculeDbRefEntry.get() or ''
      text.strip()

      if text == '':
        showWarning('Failure','Name cannot be an empty string.', parent=self)
        return

      self.align = self.molecule.findFirstAlignment()

      if not self.align:
        showWarning('Failure','No database name set for this molecule yet.', parent=self)
        return

      db = self.align.dbRef.database

      if not db:
        showWarning('Failure','No database for this alignment yet.', parent=self) # Shouldn't happen?
        return

      dbRef = db.findFirstEntry(name=text)

      if not dbRef:
        dbRef = db.newEntry(name=text)

      appData = dbRef.findFirstApplicationData(application='nmrStar',keyword='authDbAcc')

      if appData:
        for appData in dbRef.findAllApplicationData(application='nmrStar',keyword='authDbAcc'):
          dbRef.removeApplicationData(appData)

      keywds = {'application': 'nmrStar',
                'keyword':     'authDbAcc',
                'value':       'yes'}

      appData = Implementation.AppDataString(**keywds)
      dbRef.addApplicationData(appData)

      dbRef.code = text

      self.accNum = dbRef.name

      #print 'DBS: [%s] [%s]' % (db.name, self.accNum)

      self.align.__dict__['dbRef'] = dbRef
      #self.align.dbRef.delete()
      #self.align.setDbRef(dbRef)

  def setMoleculeDbRefDetails(self, event):

    if self.molecule:
      text = self.moleculeDbRefDetailsEntry.get() or ''
      text.strip()

      self.align = self.molecule.findFirstAlignment()

      if not self.align:
        showWarning('Failure','No database name set for this molecule yet.', parent=self)
        return

      db = self.align.dbRef.database

      if not db:
        showWarning('Failure','No database for this alignment yet.', parent=self) # Shouldn't happen?
        return

      if text != '':
        self.align.dbRef.details = text

  def updateMolDataAfter(self, obj=None):
    
    if self.tabbedFrame.selected == 5:
      self.after_idle(self.updateChains)
      self.after_idle(self.updateMolecules)

  def updateChains(self):

    textMatrix  = []
    objectList  = []

    if self.molSystem:
      for chain in self.molSystem.sortedChains():
        datum = [chain.code,
                 len(chain.residues),
                 chain.molecule.name,
                 chain.conformationalIsomer,
                 chain.chemExchangeState,
                 chain.physicalState]

        objectList.append(chain)
        textMatrix.append(datum)

    self.molSysTable.update(textMatrix=textMatrix,
                            objectList=objectList)

  def getMolecules(self):

    molecules = []
    if self.molSystem:
      molecules = list(set([c.molecule for c in self.molSystem.chains]))
      molecules.sort()

    return molecules

  def updateMolecules(self):

    textMatrix  = []
    objectList  = []

    molecules = self.getMolecules()

    for molecule in molecules:
      domain = mutation = ecnums = None

      molResidues = molecule.sortedMolResidues()

      seq = ''.join([mr.chemComp.code1Letter or 'X' for mr in molResidues[:25]])

      if len(molResidues) > 25:
        seq += '...'

      appData = molecule.findFirstApplicationData(application='nmrStar',keyword='Fragment')

      if appData:
        domain = appData.value

      appData = molecule.findFirstApplicationData(application='nmrStar',keyword='Mutation')

      if appData:
        mutation = appData.value

      aln = molecule.findFirstAlignment()

      dbName = accNum = dbRefdetails = None

      if aln:
        dbRef = aln.dbRef
        if dbRef:
          accNum = dbRef.name
          dbRefdetails = dbRef.details
          dbName = dbRef.database.name

      msn = molecule.findFirstMoleculeSysName(namingSystem='EC')

      if msn:
        ecnums = msn.name

      datum = [molecule.name,
               molecule.molType,
               len(molResidues),
               seq,
               molecule.isParamagnetic and 'Yes' or 'No',
               dbName,
               accNum,
               dbRefdetails,
               domain,
               ecnums,
               mutation,
               molecule.details]

      objectList.append(molecule)
      textMatrix.append(datum)

    self.moleculeTable.update(textMatrix=textMatrix,
                              objectList=objectList)

    textMatrix  = []
    objectList  = []

    chemComps = set()
    for molecule in molecules:
      for molResidue in molecule.molResidues:
        chemComp = molResidue.chemComp

        if chemComp.molType == 'other':
          chemComps.add(chemComp)
        elif not chemComp.isLinearPolymer:
          chemComps.add(chemComp)
        elif chemComp.stdChemComp and (chemComp.stdChemComp is not chemComp):
          chemComps.add(chemComp)
        elif chemComp.className == 'NonStdChemComp':
          chemComps.add(chemComp)

    chemComps = list(chemComps)
    chemComps.sort()

    for chemComp in chemComps:

      datum = [chemComp.ccpCode,
               chemComp.code3Letter,
               chemComp.molType,
               chemComp.name]

      objectList.append(chemComp)
      textMatrix.append(datum)

    self.unusualChemCompTable.update(textMatrix=textMatrix,
                                      objectList=objectList)


  # Source functions

  def selectEntryMol(self, obj, row, col):

    self.entryMolecule = obj

  def selectNatSource(self, obj, row, col):

    self.natSource = obj

    for entMol in self.entry.sortedEntryMolecules():
      if self.natSource.findFirstMolecule().findFirstEntryMolecule() == entMol:
        self.entryMolecule = entMol
        break

  def updateSourcesAfter(self, person):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateSources)

  def updateSources(self):

    if self.entry:
      if self.entryMolecule:
        if self.entryMolecule not in self.entry.entryMolecules:
          self.entryMolecule = None
      else:
        self.entryMolecule = self.entry.findFirstEntryMolecule()
    else:
      self.entryMolecule = None

    textMatrix  = []
    objectList  = []

    if self.entry:

      for entryMol in self.entry.entryMolecules:
        expSource = None
        if entryMol:
          expSource = entryMol.experimentalSource

        if expSource:
          datum = [expSource.scientificName,
                   expSource.organismName,
                   entryMol.molecule.name,
                   entryMol.productionMethod,
                   expSource.strain,
                   entryMol.vectorType,
                   expSource.plasmid,
                   expSource.cellLine,
                   expSource.variant,
                   expSource.details]

          objectList.append(entryMol)
          textMatrix.append(datum)

      self.expSourceTable.update(textMatrix=textMatrix,
                                 objectList=objectList)

      textMatrix  = []
      objectList  = []

      for entryMol in self.entry.entryMolecules:
        natSource = None
        if entryMol:
          natSource = entryMol.molecule.naturalSource

        if natSource:
          datum = [natSource.scientificName,
                   natSource.organismName,
                   entryMol.molecule.name,
                   entryMol.sourceType,
                   natSource.strain,
                   natSource.variant,
                   natSource.geneMnemonic,
                   natSource.cellLine,
                   natSource.atccNumber,
                   natSource.organ,
                   natSource.tissue,
                   natSource.cellType,
                   natSource.organelle,
                   natSource.details]

          objectList.append(natSource)
          textMatrix.append(datum)

      self.natSourceTable.update(textMatrix=textMatrix,
                                 objectList=objectList)

    else:
      self.expSourceTable.update(textMatrix=[],
                                 objectList=[])

      self.natSourceTable.update(textMatrix=[],
                                 objectList=[])
    self.waiting = False

  def addExpSource(self):

    if not self.molSystem:
      showWarning('Failure','No molSystem selected.', parent=self)
      return

    entry = self.entry

    molecules = set()
    for chain in self.molSystem.chains:
      molecule = chain.molecule
      entryMolecule = entry.findFirstEntryMolecule(molecule=molecule)
      
      if not entryMolecule:
        molecules.add(molecule)
    
    if not molecules:
      chain = self.molSystem.findFirstChain()
      
      if chain:
        molecules.add(chain.molecule)
            
    # First time round this will make at least
    #   one exp source for each molecule.

    for molecule in molecules:

      # wb104 9 Jul 2013: why this restriction??
      if len(molecule.molResidues) < 2:
        continue
       
      # wb104 9 Jul 2013: productionMethod is not part of the key
      # for EntryMolecule so the checks below are all wrong
      """
      for prodMethod0 in PRODUCTION_METHODS:
        if not entry.findFirstEntryMolecule(molecule=molecule,
                                       productionMethod=prodMethod0):
          prodMethod = prodMethod0
          break
          
      else:
        prodMethod = UNDEFINED
        while entry.findFirstEntryMolecule(molecule=molecule,
                                           productionMethod=prodMethod):
          prodMethod += '?'

      entryMolecule = entry.newEntryMolecule(molecule=molecule,
                                             productionMethod=prodMethod)
"""
      # wb104 9 Jul 2013: added below to find entryMolecule
      entryMolecule = entry.findFirstEntryMolecule(molecule=molecule)
      if not entryMolecule:
        prodMethod = PRODUCTION_METHODS[0] # arbitrary but this is a mandatory attribute so we need something
        entryMolecule = entry.newEntryMolecule(molecule=molecule,
                                               productionMethod=prodMethod)

      getExpSource(entryMolecule)

  def removeExpSource(self):

    if self.entryMolecule:
      expSource = self.entryMolecule.experimentalSource
      
      self.entryMolecule.delete()
      
      if expSource and not (expSource.entryMolecules or expSource.molecules):
        expSource.delete()      
    
      self.entryMolecule = None

    #self.waiting = False

  def setExpHostSciName(self, event):

    if self.entryMolecule:
      names = self.expHostSciNamePulldown.getObject()
      expSource = getExpSource(self.entryMolecule)
 
      if names is None:
        msg = 'Source organism scientific name:'
        sciName = askString('Input', msg, parent=self) or ''
        sciName.strip()
        
        if sciName:
           expSource.scientificName = sciName
           
           if expSource.organismName == UNDEFINED:
             expSource.organismName = sciName

      else:
        expSource.scientificName = names[0]
        
        if expSource.organismName == UNDEFINED:
          expSource.organismName = names[1] or names[0]
 
  def setExpHostName(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      name = self.expHostNameEntry.get()
      if not name:
        showError('Common Name not set','Common Name must be set', parent=self)
        return
      expSource.organismName = name

  def setExpHostMol(self, event):

    entryMolecule = self.entryMolecule
    if entryMolecule:
      molecule = self.expHostMolPulldown.getObject() or None
      
      if molecule:
        prodMethod = entryMolecule.productionMethod
        sourceType = entryMolecule.sourceType
        vectorType = entryMolecule.vectorType
        expSource  = getExpSource(entryMolecule)
        
        entryMolecule.delete()
        entryMolecule = self.entry.newEntryMolecule(molecule=molecule,
                                                    experimentalSource=expSource,
                                                    productionMethod=prodMethod,
                                                    vectorType=vectorType,
                                                    sourceType=sourceType)

        self.entryMolecule = None

  def setExpHostProdMethod(self, event):

    if self.entryMolecule:
      prodMethod = self.expHostProdMethodPulldown.getObject()

      if prodMethod is None:
        msg = 'Non standard production method:'
        prodMethod = askString('Input', msg, parent=self) or ''
        prodMethod.strip()

      self.entryMolecule.productionMethod = prodMethod

  def setExpHostStrain(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      expSource.strain = self.expHostStrainEntry.get() or None

  def setExpHostVecType(self, event):

    if self.entryMolecule:
      vecType = self.expHostVecTypePulldown.getObject()

      if vecType is None:
        msg = 'Non standard vector type:'
        vecType = askString('Input', msg, parent=self) or ''
        vecType.strip()

      self.entryMolecule.vectorType = vecType or None

  def setExpHostVecName(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      expSource.plasmid = self.expHostVecNameEntry.get() or None

  def setExpHostCellLine(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      expSource.cellLine = self.expHostCellLineEntry.get() or None

  def setExpHostVariant(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      expSource.variant = self.expHostVariantEntry.get() or None

  def setExpHostDetails(self, event):

    if self.entryMolecule:
      expSource = getExpSource(self.entryMolecule)
      expSource.details = self.expHostDetailsEntry.get() or None

  def getExpHostSciName(self, entryMolecule):

    index = 0
    names = []
    organisms = []
    #categories = []

    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      sciName0 = expSource.scientificName
      
      for sciName in EXP_ORGANISMS:
          names.append(sciName)
          organisms.append((sciName, sciName))
          
          if sciName == sciName0:
            index = len(organisms)-1

      names.append('<Other>')
      organisms.append(None)
      #categories.append(None)
      
    self.expHostSciNamePulldown.setup(names, organisms, index)

  def getExpHostName(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.organismName

    self.expHostNameEntry.set(value)

  def getExpHostMol(self, entryMolecule):

    molecules = self.getMolecules()

    index = 0
    if self.entryMolecule:
      index = molecules.index(self.entryMolecule.molecule)

    names = [ mol.name for mol in molecules ]

    self.expHostMolPulldown.setup(names, molecules, index)

  def getExpHostProdMethod(self, entryMolecule):

    if not entryMolecule:
      return

    if not hasattr(self, 'prodMethods'):
      self.prodMethods = PRODUCTION_METHODS
      self.prodMethods.sort()

    cProdMethod = entryMolecule.productionMethod

    if cProdMethod not in self.prodMethods:
      self.prodMethods.append(cProdMethod)
      self.prodMethods.sort()

    prodMethods = self.prodMethods[:]

    # Filter out already used methods
    # This didn't seem to work for me (CJP)???

    #entry = entryMolecule.entry
    #for prodMethod in self.prodMethods:
    #  if not entry.findFirstEntryMolecule(molecule=entryMolecule.molecule,
    #                                      productionMethod=prodMethod):
    #    prodMethods.append(prodMethod)

    #print "METH: [%s] [%s]" % (prodMethods, cProdMethod)
    
    index = prodMethods.index(cProdMethod)
    names = prodMethods[:] + ['<Other>',]
    prodMethods.append(None)

    self.expHostProdMethodPulldown.setup(names, prodMethods, index)

  def getExpHostStrain(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.strain

    self.expHostStrainEntry.set(value)

  def getExpHostVecType(self, expSource):

    if not hasattr(self, 'vecTypes'):
      self.vecTypes = VECTOR_TYPES
      self.vecTypes.sort()
      self.vecTypes.append(None)

    cVecType = None
    if self.entryMolecule:
      cVecType = self.entryMolecule.vectorType

    if cVecType:
      if cVecType not in self.vecTypes:
        self.vecTypes.insert(0, cVecType)
        self.vecTypes.remove(None)
        self.vecTypes.sort()
        self.vecTypes.append(None)

    index = self.vecTypes.index(cVecType)

    names = self.vecTypes[:-1] + ['<Other>',]

    self.expHostVecTypePulldown.setup(names, self.vecTypes, index)

  def getExpHostVecName(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.plasmid

    self.expHostVecNameEntry.set(value)

  def getExpHostCellLine(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.cellLine

    self.expHostCellLineEntry.set(value)

  def getExpHostVariant(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.variant

    self.expHostVariantEntry.set(value)

  def getExpHostDetails(self, entryMolecule):

    value = None
    if entryMolecule:
      expSource = getExpSource(entryMolecule)
      value = expSource.details

    self.expHostDetailsEntry.set(value)

  def addNatSources(self):

    if not self.molSystem:
      showWarning('Failure','No molSystem selected.', parent=self)
      return

    entry = self.entry

    molecules = set()
    for chain in self.molSystem.chains:
      molecule = chain.molecule
      
      if not molecule.naturalSource:
        molecules.add(molecule)
                
    # First time round this will make at least
    # one natural source for each molecule
    
    tax = getTaxonomy(entry.root)
    for molecule in molecules:

      if len(molecule.molResidues) < 2:
        continue

      #for natSource in tax.naturalSources:
      #  if natSource.molecules:
      #    molecule.naturalSource = natSource
      #    break
      #    
      #else:
      
      natSource = tax.newNaturalSource(organismName=UNDEFINED,
                                       scientificName=UNDEFINED)
      molecule.naturalSource = natSource

  def removeNatSource(self):

    if self.natSource:
      self.natSource.delete()
      self.natSource = None

    #self.waiting = False

  def setNatSourceSciName(self, event):

    if self.natSource:
      names = self.natSourceSciNamePulldown.getObject()

      if names and names[0]:

        for sciName, commonName in ORGANISMS['Viruses']:
          if names[0] == sciName:
            self.entryMolecule.sourceType = 'virus'

      if names is None:
        msg = 'Source organism scientific name:'
        sciName = askString('Input', msg, parent=self) or ''
        sciName.strip()
        
        if sciName:
           self.natSource.scientificName = sciName
           
           if self.natSource.organismName == UNDEFINED:
             self.natSource.organismName = sciName

      else:
        self.natSource.scientificName = names[0]
        
        if self.natSource.organismName == UNDEFINED:
          self.natSource.organismName = names[1] or names[0]

  def setNatSourceName(self, event):

    if self.natSource:
      name = self.natSourceNameEntry.get()
      if not name:
        showError('Common Name not set','Common Name must be set', parent=self)
        return
      self.natSource.organismName = self.natSourceNameEntry.get()

  def setNatSourceMol(self, event):

    if self.natSource:
      molecule = self.natSourceMolPulldown.getObject() or None
      if molecule:
        self.natSource.molecule = molecule

  def setNatSourceOrgType(self, event):

    if self.natSource:
      orgType = self.natSourceOrgTypePulldown.getObject()

      if orgType is None:
        msg = 'Non standard organism type:'
        orgType = askString('Input', msg, parent=self) or ''
        orgType.strip()

      self.entryMolecule.sourceType = orgType or None

  def setNatSourceStrain(self, event):

    if self.natSource:
      self.natSource.strain = self.natSourceStrainEntry.get() or None

  def setNatSourceVariant(self, event):

    if self.natSource:
      self.natSource.variant = self.natSourceVariantEntry.get() or None

  def setNatSourceGeneName(self, event):

    if self.natSource:
      self.natSource.geneMnemonic = self.natSourceGeneNameEntry.get() or None

  def setNatSourceCellLine(self, event):

    if self.natSource:
      self.natSource.cellLine = self.natSourceCellLineEntry.get() or None

  def setNatSourceATCCNum(self, event):

    if self.natSource:
      self.natSource.atccNumber = self.natSourceATCCNumEntry.get() or None

  def setNatSourceOrgan(self, event):

    if self.natSource:
      self.natSource.organ = self.natSourceOrganEntry.get() or None

  def setNatSourceTissue(self, event):

    if self.natSource:
      self.natSource.tissue = self.natSourceTissueEntry.get() or None

  def setNatSourceCellType(self, event):

    if self.natSource:
      self.natSource.cellType = self.natSourceCellTypeEntry.get() or None

  def setNatSourceOrganelle(self, event):

    if self.natSource:
      self.natSource.organelle = self.natSourceOrganelleEntry.get() or None

  def setNatSourceDetails(self, event):

    if self.natSource:
      self.natSource.details = self.natSourceDetailsEntry.get() or None

  def getNatSourceSciName(self, natSource):

    index = 0
    names = []
    organisms = []
    categories = []

    if self.natSource:
      sciName0 = self.natSource.scientificName
      
      for kingdom in sorted(ORGANISMS.keys() ):
        for sciName, commonName in ORGANISMS[kingdom]:
          names.append(sciName)
          organisms.append((sciName, commonName))
          categories.append(kingdom)
          if sciName == sciName0:
            index = len(organisms)-1

      names.append('<Other>')
      organisms.append(None)
      categories.append(None)
      
    self.natSourceSciNamePulldown.setup(names, organisms, index, categories=categories)

  def getNatSourceName(self, natSource):

    value = None
    if natSource:
      value = natSource.organismName

    self.natSourceNameEntry.set(value)

  def getNatSourceMol(self, natSource):

    molecules = self.getMolecules()

    index = 0
    if self.entryMolecule:
      index = molecules.index(self.entryMolecule.molecule)

    names = [ mol.name for mol in molecules ]

    self.natSourceMolPulldown.setup(names, molecules, index)

  def getNatSourceOrgType(self, natSource):

    if not hasattr(self, 'orgTypes'):
      self.orgTypes = ORGANISM_TYPES
      self.orgTypes.sort()
      self.orgTypes.append(None)

    cOrgType = None
    if self.entryMolecule:
      cOrgType = self.entryMolecule.sourceType

    if cOrgType:
      if cOrgType not in self.orgTypes:
        self.orgTypes.insert(0, cOrgType)
        self.orgTypes.remove(None)
        self.orgTypes.sort()
        self.orgTypes.append(None)

    index = self.orgTypes.index(cOrgType)

    names = self.orgTypes[:-1] + ['<Other>',]

    self.natSourceOrgTypePulldown.setup(names, self.orgTypes, index)

  def getNatSourceStrain(self, natSource):

    value = None
    if natSource:
      value = natSource.strain

    self.natSourceStrainEntry.set(value)

  def getNatSourceVariant(self, natSource):

    value = None
    if natSource:
      value = natSource.variant

    self.natSourceVariantEntry.set(value)

  def getNatSourceGeneName(self, natSource):

    value = None
    if natSource:
      value = natSource.geneMnemonic

    self.natSourceGeneNameEntry.set(value)

  def getNatSourceCellLine(self, natSource):

    value = None
    if natSource:
      value = natSource.cellLine

    self.natSourceCellLineEntry.set(value)

  def getNatSourceATCCNum(self, natSource):

    value = None
    if natSource:
      value = natSource.atccNumber

    self.natSourceATCCNumEntry.set(value)

  def getNatSourceOrgan(self, natSource):

    value = None
    if natSource:
      value = natSource.organ

    self.natSourceOrganEntry.set(value)

  def getNatSourceTissue(self, natSource):

    value = None
    if natSource:
      value = natSource.tissue

    self.natSourceTissueEntry.set(value)

  def getNatSourceCellType(self, natSource):

    value = None
    if natSource:
      value = natSource.cellType

    self.natSourceCellTypeEntry.set(value)

  def getNatSourceOrganelle(self, natSource):

    value = None
    if natSource:
      value = natSource.organelle

    self.natSourceOrganelleEntry.set(value)

  def getNatSourceDetails(self, natSource):

    value = None
    if natSource:
      value = natSource.details

    self.natSourceDetailsEntry.set(value)
    

  # Experiment functions

  def checkExperimentRemoval(self, experiment, verbose=False):

    for spectrum in experiment.dataSources:
      for peakList in spectrum.peakLists:
        if peakList in self.entry.peakLists:
          if verbose:
            msg = "Cannot toggle experiment off while it's "
            msg += "peak lists are selected for deposition."
            showWarning('Warning', msg, parent=self)
          return False

    # Note (also see checkEntryConsistency sub-routine):
    #   this will tell you if the shiftList for the experiment is in the entry and so can't be toggled off.

    #if experiment.shiftList:
    #  if experiment.shiftList in self.entry.findAllMeasurementLists(className='ShiftList'):
    #    if verbose:
    #      msg = "Cannot toggle experiment off while it's "
    #      msg += "shift list is selected for deposition."
    #      showWarning('Warning', msg, parent=self)
    #    return False

    constraintSets = set()
    for structureAnalysis in self.entry.structureAnalyses:
      rSet = structureAnalysis.nmrConstraintStore
      if rSet:
        constraintSets.add(rSet)

    for structureAnalysis in self.entry.structureGenerations:
      rSet = structureAnalysis.nmrConstraintStore
      if rSet:
        constraintSets.add(rSet)

    for rSet in constraintSets:
      for constraintList in rSet.constraintLists:
        if experiment in constraintList.experiments:
          if verbose:
            msg = "Cannot toggle experiment off as it is associated with "
            msg += "the structural restraints selected for deposition."
            showWarning('Warning', msg, parent=self)
          return False

    return True

  def toggleExperiment(self, experiment):

    if experiment in self.entry.experiments:
      present = True
    else:
      present = False

    if present:
      if not self.checkExperimentRemoval(experiment, verbose=True):
        return

      self.entry.removeExperiment(experiment)

    else:
      self.entry.addExperiment(experiment)

  def selectExperiment(self, experiment, row, col):

    self.experiment = experiment

  def addExperiment(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    msg = 'Name of experiment'
    name = askString('Query', msg, parent=self) or ''
    name.strip()

    #msg = 'No of dimensions'
    #numDim = askInteger('Query', msg, parent=self) or 2
    numDim = 2

    if not self.nmrProject:
      self.nmrProject = getNmrProject(self.project)

    self.experiment = self.nmrProject.newExperiment(name=name, numDim=numDim)

    self.entry.addExperiment(self.experiment)

  def removeExperiment(self):

    if self.experiment:
      self.experiment.delete()
      self.experiment = None

  def setExperimentName(self, event):

    if self.experiment:
      self.experiment.name = self.experimentNameEntry.get() or self.experiment.name

  def setExperimentState(self, event):

    if self.experiment:
      state = self.experimentStatePulldown.getObject()

      if state is None:
        msg = 'Other experiment state:'
        state = askString('Input', msg, parent=self) or ''
        state.strip()

      self.experiment.sampleState = state or None

  def setExperimentSample(self, event):

    if self.experiment:

      sample = self.experimentSamplePulldown.getObject()

      if sample:
        self.experiment.setSample(sample)
      else:
        self.experiment.sample = None

  def setExperimentSampleCond(self, event):

    if self.experiment:

      scs = self.experimentSampleCondPulldown.getObject()

      if scs:
        self.experiment.setSampleConditionSet(scs)
      else:
        self.experiment.sampleConditionSet = None

  def setExperimentSpectrometer(self, event):

    if self.experiment:

      spectrometer = self.experimentSpecPulldown.getObject()

      if spectrometer:
        self.experiment.setSpectrometer(spectrometer)
      else:
        self.experiment.spectrometer = None

  def setExperimentProbe(self, event):

    if self.experiment:

      probe = self.experimentProbePulldown.getObject()

      if probe:
        self.experiment.setProbe(probe)
      else:
        self.experiment.probe = None

  def setExpChemShiftList(self, event):

    if self.experiment:

      shiftList = self.expChemShiftListPulldown.getObject()

      if shiftList:
        self.experiment.setShiftList(shiftList)

      else:
        self.experiment.shiftList = None

  def setExpShiftRefs(self, obj):
  
    if self.experiment:
      if obj is None:
        self.experimentTable.keyPressEscape()
      else:
        shiftRefs = self.experiment.nmrProject.sortedShiftReferences()
        values = self.expShiftRefSelect.get()
        selectedRefs = [shiftRefs[i] for i in range(len(values)) if values[i]]
        self.experiment.setShiftReferences(selectedRefs)
        self.experimentTable.keyPressEscape()
  
  def getExpShiftRefs(self, experiment):
  
    shiftRefs = self.nmrProject.sortedShiftReferences()
    names  = []
    values = []
    for shiftReference in shiftRefs:
      data = (shiftReference.serial,
              shiftReference.isotopeCode,
              shiftReference.molName)
      names.append('%d:%s:%s' % data)
      
      if shiftReference in experiment.shiftReferences:
        values.append(True)
      else:
        values.append(False)  
        
    self.expShiftRefSelect.set(values=values,options=names)

  def getExperimentName(self, experiment):

    text = ''
    if experiment:
      text = experiment.name

    self.experimentNameEntry.set(text)

  def getExperimentState(self, experiment):

    if not hasattr(self, 'expStates'):
      self.expStates = EXP_STATES
      #self.expStates.sort()
      self.expStates.append(None)

    cExpState = None

    if experiment:
      cExpState = experiment.sampleState

    index = 0

    if cExpState:
      if cExpState not in self.expStates:
        self.expStates.insert(0, cExpState)
        #self.expStates.remove(None)
        #self.expStates.sort()
        #self.expStates.append(None)

      index = self.expStates.index(cExpState)

    names = self.expStates[:-1] + ['<Other>',]

    self.experimentStatePulldown.setup(names, self.expStates, index)

  def getExperimentSample(self, experiment):

    sampleStore = getSampleStore(self.project)

    index = 0
    samples = []
    names = []
    for sample in sampleStore.sortedAbstractSamples():
      names.append(sample.name)
      samples.append(sample)

    names.append('<None>')
    samples.append(None)

    if experiment.sample:
      index = samples.index(experiment.sample)

    self.experimentSamplePulldown.setup(names, samples, index)

  def getExperimentSampleCond(self, experiment):

    if not self.nmrProject:
      self.nmrProject = getNmrProject(self.project)

    index = 0
    scsList = []
    names = []
    for scs in self.nmrProject.sortedSampleConditionSets():
      names.append(scs.name)
      scsList.append(scs)

    names.append('<None>')
    scsList.append(None)

    if experiment.sampleConditionSet:
      index = scsList.index(experiment.sampleConditionSet)

    self.experimentSampleCondPulldown.setup(names, scsList, index)

  def getExperimentSpectrometer(self, experiment):

    instrumentStore = getInstrumentStore(self.project)

    if not instrumentStore:
      return

    getInstruments = instrumentStore.findAllInstruments

    spectrometers = list(getInstruments(className='NmrSpectrometer') )
    spectrometers.sort()
    
    index = 0
    names = []
    for spectrometer in spectrometers:
      names.append(spectrometer.name)

    names.append('<None>')
    spectrometers.append(None)

    if experiment.spectrometer:
      index = spectrometers.index(experiment.spectrometer)

    self.experimentSpecPulldown.setup(names, spectrometers, index)

  def getExperimentProbe(self, experiment):

    instrumentStore = getInstrumentStore(self.project)
    getInstruments = instrumentStore.findAllInstruments

    probes = list(getInstruments(className='NmrProbe') )
    probes.sort()
    
    index = 0
    names = []
    for probe in probes:
      names.append(probe.name)

    names.append('<None>')
    probes.append(None)

    if experiment.probe:
      index = probes.index(experiment.probe)

    self.experimentProbePulldown.setup(names, probes, index)

  def getExpChemShiftList(self, experiment):

    if not self.entry:
      return

    shiftLists = list(self.entry.findAllMeasurementLists(className='ShiftList') )

    if not shiftLists:
      showWarning('Failure','You need to select a shift list for deposition in the "NMR Data" tab.', parent=self)
      return

    shiftLists.sort()
    
    index = 0
    names = []
    for sl in shiftLists:
      names.append(str(sl.serial) + ':' + sl.name)

    names.append('<None>')
    shiftLists.append(None)

    if experiment and experiment.shiftList and shiftLists and experiment.shiftList in shiftLists:
      index = shiftLists.index(experiment.shiftList)

    self.expChemShiftListPulldown.setup(names, shiftLists, index)

  def submitAllExperiments(self):

    if self.entry and self.nmrProject:
      for experiment in self.nmrProject.sortedExperiments():
        if experiment not in self.entry.experiments:
          self.entry.addExperiment(experiment)

  def submitMinExperiments(self):   

    if self.entry:
      for experiment in self.entry.experiments:
        if self.checkExperimentRemoval(experiment):
          self.entry.removeExperiment(experiment)

  def updateExperimentsAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateExperiments)

  def getStarExperimentName(self,experiment):

    bmrbName = None

    if experiment.refExperiment:
      refExperiment = experiment.refExperiment

      bmrbSysName = refExperiment.findFirstSystematicName(namingSystem = 'BMRB')

      if bmrbSysName:
        bmrbName = bmrbSysName.name

      if not bmrbName:
        bmrbName = refExperiment.name

        if not bmrbName and refExperiment.nmrExpPrototype:
          bmrbName = refExperiment.nmrExpPrototype.name

          if not bmrbName:
            bmrbName = refExperiment.nmrExpPrototype.synonym

      if bmrbName and not (bmrbSysName and bmrbSysName.name):
        bmrbName = experiment.name + ' (' + bmrbName + ')'

    if not bmrbName:
      bmrbName = experiment.name

    return bmrbName

  def updateExperiments(self):

    self.updateSampleConditionSets()

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if self.entry and self.nmrProject:

      nCols = len(self.experimentTable.headingList)

      for experiment in self.nmrProject.sortedExperiments():

        if experiment.spectrometer:
          spectrometer = experiment.spectrometer.name
        else:
          spectrometer = None

        if experiment.probe:
          probe = experiment.probe.name
        else:
          probe = None

        if experiment.sample:
          sample = experiment.sample.name
        else:
          sample = None

        if experiment.sampleConditionSet:
          scs = experiment.sampleConditionSet.name
        else:
          scs = None

        if experiment.shiftList:
          shiftList = str(experiment.shiftList.serial) + ':' + experiment.shiftList.name
        else:
          shiftList = None

        if experiment in self.entry.experiments:
          submit = 'Yes'
          colors = [NICE_GREEN] * nCols
        else:
          submit = 'No'
          colors = [None] * nCols

        shiftRefs = []
        for shiftReference in experiment.sortedShiftReferences():
          data = (shiftReference.isotopeCode,
                  shiftReference.molName)
          shiftRefs.append('%s:%s' % data)
 
        shiftRef = ','.join(shiftRefs)

        bmrbName = self.getStarExperimentName(experiment)

        datum = [experiment.name,
                 submit,
                 bmrbName,
                 experiment.sampleState,
                 sample,
                 scs,
                 shiftRef,
                 experiment.rawData and 'Yes' or 'No',
                 spectrometer,
                 probe,
                 shiftList]

        objectList.append(experiment)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    self.experimentTable.update(textMatrix=textMatrix,
                                objectList=objectList,
                                colorMatrix=colorMatrix)

    self.waiting = False

  """
  Experiment set attrs:

  sampleState
  sampleVolume
  volumeUnit
  sample
  """


  # Sample functions

  def selectSample(self, sample, row, col):

    self.sample = sample

    self.updateSampleComponentsAfter()

  def addSample(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    msg = 'Name of sample'
    name = askString('Query', msg, parent=self) or ''
    name.strip()

    if name == '':
      showWarning('Failure','Name cannot be an empty string.', parent=self)
      return

    if len(name.split()) > 1:
      showWarning('Failure','Name cannot contain whitespace.', parent=self)
      return

    sampleStore = getSampleStore(self.project)

    if sampleStore.findFirstAbstractSample(name=name):
      showWarning('Failure','Name already in use.', parent=self)
      return

    classification = getClassification(self.project)

    if not hasattr(self, 'sampleCategory'):
      self.sampleCategory = classification.findFirstSampleCategory(name='deposition')
      if not self.sampleCategory:
        self.sampleCategory = classification.newSampleCategory(name='deposition')

    sampleStore.newSample(name=name, sampleCategories=[self.sampleCategory])

  def removeSample(self):  

    if self.sample:
      self.sample.delete()
      self.sample = None

  """
  WV 22/07/2010 - removed, should be handled by CCPN to NMR-STAR export
  def setSampleState(self, event):

    if self.sample:
      appData = self.sample.findFirstApplicationData(application='nmrStar',keyword='sampleState')

      if appData:
        for appData in self.sample.findAllApplicationData(application='nmrStar',keyword='sampleState'):
          self.sample.removeApplicationData(appData)

      state = self.sampleStatePulldown.getObject()

      if state is None:
        msg = 'Other sample state:'
        state = askString('Input', msg, parent=self) or ''
        state.strip()

      if state:
        keywds = {'application': 'nmrStar',
                  'keyword':     'sampleState',
                  'value':       state}

        appData = Implementation.AppDataString(**keywds)
        self.sample.addApplicationData(appData)
  """

  def setSampleAmount(self, event):

    if self.sample:
      self.sample.initialAmount = self.sampleAmountEntry.get() or None

  def setSampleUnit(self, event):

    if self.sample:
      unit = self.sampleUnitPulldown.getObject()
      self.sample.amountUnit = unit

  def setSampleIonic(self, event):

    if self.sample:
      ionicStr = self.sampleIonicEntry.get() or None
      self.sample.ionicStrength = ionicStr

  def setSamplePh(self, event):

    if self.sample:
      ph = self.samplePhEntry.get() or None
      self.sample.ph = ph

  def setSampleSolvent(self, event):

    if self.sample:
      appData = self.sample.findFirstApplicationData(application='nmrStar',keyword='solventSys')

      if appData:
        for appData in self.sample.findAllApplicationData(application='nmrStar',keyword='solventSys'):
          self.sample.removeApplicationData(appData)

      solventSys = self.sampleSolventPulldown.getObject()

      if solventSys is None:
        msg = 'Other solvent syatem:'
        solventSys = askString('Input', msg, parent=self) or ''
        solventSys.strip()

      if solventSys:
        keywds = {'application': 'nmrStar',
                  'keyword':     'solventSys',
                  'value':       solventSys}

        appData = Implementation.AppDataString(**keywds)
        self.sample.addApplicationData(appData)

  def setSampleDetails(self, event):

    if self.sample:
      details = self.sampleDetailsEntry.get() or ''
      if details != '':
        self.sample.details = details

  """
  def getSampleState(self, sample):

    if not hasattr(self, 'sampStates'):
      self.sampStates = SAMPLE_STATES
      #self.sampStates.sort()
      self.sampStates.append(None)

    value = None

    if sample:
      value = None

      appData = sample.findFirstApplicationData(application='nmrStar',keyword='sampleState')

      if appData:
        value = appData.value

    index = 0

    if value:
      if value not in self.sampStates:
        self.sampStates.insert(0, value)
        #self.sampStates.remove(None)
        #self.sampStates.sort()
        #self.sampStates.append(None)

      index = self.sampStates.index(value)

    names = self.sampStates[:-1] + ['<Other>',]

    self.sampleStatePulldown.setup(names, self.sampStates, index)
  """

  def getSampleAmount(self, sample):

    value = None
    if sample:
      value = sample.initialAmount

    self.sampleAmountEntry.set(value)

  def getSampleUnit(self, sample):

    unit = sample.amountUnit

    if not unit:
      return

    names = SAMPLE_UNITS

    index = names.index(unit)
    self.sampleUnitPulldown.setup(names, names, index)

  def getSampleIonic(self, sample):

    value = None
    if sample:
      value = sample.ionicStrength

    self.sampleIonicEntry.set(value)

  def getSamplePh(self, sample):

    value = None
    if sample:
      value = sample.ph

    self.samplePhEntry.set(value)

  def getSampleSolvent(self, sample):

    if not hasattr(self, 'solventSys'):
      self.solventSys = SOLVENT_SYSTEMS
      self.solventSys.sort()
      self.solventSys.append(None)

    value = None

    appData = sample.findFirstApplicationData(application='nmrStar',keyword='solventSys')

    if appData:
      value = appData.value

    index = 0

    if value:
      if value not in self.solventSys:
        self.solventSys.insert(0, value)
        self.solventSys.remove(None)
        self.solventSys.sort()
        self.solventSys.append(None)

      index = self.solventSys.index(value)

    names = self.solventSys[:-1] + ['<Other>',]

    self.sampleSolventPulldown.setup(names, self.solventSys, index)

  def getSampleDetails(self, sample):

    text = ''
    if sample:
      text = sample.details

    self.sampleDetailsEntry.set(text)

  def updateSamplesAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateSamples)

  def updateSamples(self, obj=None):

    textMatrix  = []
    objectList  = []

    if self.project:
      sampleStore = getSampleStore(self.project)

      for i, sample in enumerate(sampleStore.sortedAbstractSamples()):

        if sample.className != 'Sample':
          continue

        sampleState = None

        appData = sample.findFirstApplicationData(application='nmrStar',keyword='sampleState')
        if appData:
          sampleState = appData.value

        solventSys = None

        appData = sample.findFirstApplicationData(application='nmrStar',keyword='solventSys')
        if appData:
          solventSys = appData.value

        datum = [i+1,
                 sample.name,
                 ##sampleState,
                 #sample.initialAmount,
                 #sample.amountUnit,
                 #sample.ionicStrength,
                 #sample.ph,
                 solventSys,
                 sample.details]

        objectList.append(sample)
        textMatrix.append(datum)

    self.sampleTable.update(textMatrix=textMatrix,
                            objectList=objectList)

    self.waiting = False


  # Sample component functions

  def selectSampleComponent(self, obj, row, col):

    self.sampleComponent = obj

  def addPolySampleComponents(self):

    if not self.sample:
      showWarning('Failure','No sample selected in top frame.', parent=self)
      return

    if not self.molSystem:
      showWarning('Failure','No molSystem selected.', parent=self)
      return

    refSampleCompStore = getRefSampleComponentStore(self.project)

    molecules = self.getMolecules()

    for mol in molecules:
      if not mol.isFinalised:
        mol.isFinalised = True

      sampCompFlag = False

      for sampComp in self.sample.sortedSampleComponents():
        if sampComp.refComponent.className == 'MolComponent':
          if sampComp.refComponent.molecule == mol:
            sampCompFlag = True

      if sampCompFlag:
        continue

      """
      # wb104 8 Jul 2013: if the same sample component is used then editing the
      # Isotopic Labelling of one, for example, means the same happens to all
      # the others, which is presumably normally not intended.

      molComp = refSampleCompStore.findFirstComponent(className='MolComponent',
                                                      name=mol.name + ':' + self.sample.name)
"""
      molComp = None  # wb104: alternative to above

      if not molComp:
        # wb104 9 Jul 2013: make sure key is unique
        name0 = mol.name + ':' + self.sample.name
        name = name0 = name0[:80] # cannot be more than 80 chars
        nn = 1
        while refSampleCompStore.findFirstComponent(className='MolComponent', name=name):
          nn += 1
          name = name0[:72] + " " + nn
        molComp = refSampleCompStore.newMolComponent(name=name,
                                                     details=mol.name,
                                                     molecule=mol)

      sampComp = self.sample.newSampleComponent(refComponent=molComp)

  def addOtherSampleComponent(self):

    if not self.sample:
      showWarning('Failure','No sample selected in top frame.', parent=self)
      return

    msg = 'Other sample component'
    sampCompMol = askString('Input', msg, parent=self) or ''
    sampCompMol.strip()

    if sampCompMol == '':
      showWarning('Failure','Name cannot be an empty string.', parent=self)
      return

    refSampleCompStore = getRefSampleComponentStore(self.project)

    for comp in refSampleCompStore.sortedComponents():
      if sampCompMol == comp.name:
        # TODO: warning
        return

    molComp = refSampleCompStore.newMolComponent(name=sampCompMol + ':' + self.sample.name, # TODO: better?
                                                 details = sampCompMol)

    sampComp = self.sample.newSampleComponent(refComponent=molComp)

  def removeSampleComponent(self):

    if self.sampleComponent:

      msg = 'Really remove sample component?'

      if showOkCancel('query',msg):
        self.sampleComponent.delete()
        self.sampleComponent = None

  def updateSampleComponentsAfter(self, *opt):

    #if self.waiting:
    #  return
    #else:
    #  self.waiting = True

    self.after_idle(self.updateSampleComponents)

  def updateSampleComponents(self, obj=None):

    from ccpnmr.eci.IsotopeLabeling import getStarIsotopeLabeling

    textMatrix  = []
    objectList  = []

    if self.project and self.sample:

      for sampComp in self.sample.sortedSampleComponents():

        molName = ''

        if sampComp.refComponent.className == 'MolComponent':
          if sampComp.refComponent.molecule:
            molName = sampComp.refComponent.molecule.name
          else:
            molName = sampComp.refComponent.details
        else:
          molName = sampComp.refComponent.details

        bmrbLabel = getStarIsotopeLabeling(sampComp.refComponent)

        datum = [molName,
                 bmrbLabel,
                 sampComp.concentration,
                 sampComp.concentrationError,
                 sampComp.concentrationUnit]

        objectList.append(sampComp)
        textMatrix.append(datum)

    self.sampleComponentTable.update(textMatrix=textMatrix,
                                     objectList=objectList)

    #self.waiting = False

  def setSampleCompMol(self, event):

    if not self.sample:
      return # TODO: warning here

    if self.sampleComponent:
      sampCompMol = self.sampleCompMolPulldown.getObject()

      if sampCompMol is None:
        msg = 'Other sample component'
        sampCompMol = askString('Input', msg, parent=self) or ''
        sampCompMol.strip()

      molecules = self.getMolecules()

      for mol in molecules:
        if sampCompMol == mol.name:
          if not mol.isFinalised:
            mol.isFinalised = True

          if self.sampleComponent.refComponent.className == 'MolComponent':
            self.sampleComponent.refComponent.setMolecule(mol)
            break
      else:
        #self.sampleComponent.refComponent.name = sampCompMol + self.sample.name # TODO: better
        self.sampleComponent.refComponent.details = sampCompMol

  def setSampleCompIsotope(self, event):

    from ccpnmr.eci.IsotopeLabeling import getStarIsotopeLabeling, makeLabelObjects

    if not self.sample:
      return # TODO: warning here

    if self.sampleComponent:
      refComp = self.sampleComponent.refComponent

      currentBmrbLabelName = getStarIsotopeLabeling(refComp)

      bmrbLabelName = self.sampleCompIsotopePulldown.getObject()

      if bmrbLabelName is None:
        msg = 'Other isotope labeling scheme'
        bmrbLabelName = askString('Input', msg, parent=self) or ''
        bmrbLabelName.strip()

      #print 'LABEL: [%s]' % sampCompIsotope

      if bmrbLabelName and bmrbLabelName != currentBmrbLabelName:

        makeLabelObjects(self.project, refComp, bmrbLabelName)

        self.updateSampleComponentsAfter()

  def getSampleCompIsotope(self, sampComp):

    from ccpnmr.eci.IsotopeLabeling import getStarIsotopeLabeling

    cLabel = getStarIsotopeLabeling(sampComp.refComponent)

    names = []
    labels = []
    categories = []
    index = 0

    for label in ISOTOPE_LABEL_LIST:
      names.append(label[1])
      labels.append(label[1])
      categories.append(label[0])
          
      if cLabel == label[1]:
        index = len(labels)-1

    names.append('<Other>')
    labels.append(None)
    categories.append(None)
      
    self.sampleCompIsotopePulldown.setup(names, labels, index, categories=categories)

  def setSampleCompConc(self, event):

    if self.sampleComponent:
      conc = self.sampleCompConcEntry.get() or None
      self.sampleComponent.concentration = conc

  def setSampleCompConcErr(self, event):

    if self.sampleComponent:
      err = self.sampleCompConcErrEntry.get() or None
      self.sampleComponent.concentrationError = err

  def setSampleCompConcUnit(self, event):

    if self.sampleComponent:
      unit = self.sampleCompConcUnitPulldown.getObject()
      self.sampleComponent.concentrationUnit = unit

  def getSampleCompMol(self, sampleComponent):

    if not self.sample or self.sampleComponent:
      return

    molecules = self.getMolecules()

    molNames = [ mol.name for mol in molecules ]

    sampCompMols = molNames
    sampCompMols.sort()

    for sampComp in self.sample.sortedSampleComponents():
      if sampComp.refComponent.className == 'MolComponent':
        if sampComp.refComponent.molecule is not None:
          continue

      sampCompMols.append(sampComp.refComponent.details)

    sampCompMols.append(None)

    cMolName = None
    if sampleComponent:
      if sampleComponent.refComponent.className == 'MolComponent':
        if sampleComponent.refComponent.molecule:
          cMolName = sampleComponent.refComponent.molecule.name
        else:
          cMolName = sampleComponent.refComponent.details
      else:
        cMolName = sampleComponent.refComponent.details


    if cMolName:
      if cMolName not in sampCompMols:
        sampCompMols.insert(0, cMolName)

    index = sampCompMols.index(cMolName)

    names = sampCompMols[:-1] + ['<Other>',]

    self.sampleCompMolPulldown.setup(names, sampCompMols, index)

  def getSampleCompConc(self, sampComp):

    value = None
    if sampComp:
      value = sampComp.concentration

    self.sampleCompConcEntry.set(value)

  def getSampleCompConcErr(self, sampComp):

    value = None
    if sampComp:
      value = sampComp.concentrationError

    self.sampleCompConcErrEntry.set(value)

  def getSampleCompConcUnit(self, sampComp):

    unit = sampComp.concentrationUnit

    if not unit:
      return

    names = SAMPLE_COMP_UNITS

    index = names.index(unit)
    self.sampleCompConcUnitPulldown.setup(names, names, index)


  # Spectrometer functions

  def selectSpectrometer(self, obj, row, col):

    self.spectrometer = obj

  def addSpectrometer(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    instrumentStore = getInstrumentStore(self.project)

    name = ''
    while not name:
      name = askString('Text Entry','Enter NMR Spectrometer Name')

    spectrometer = instrumentStore.newNmrSpectrometer(name=name)

  def removeSpectrometer(self):

    if self.spectrometer:

      msg = 'Really remove specification of NMR Spectrometer %s?' %  self.spectrometer.name

      if showOkCancel('query',msg):
        self.spectrometer.delete()
        self.spectrometer = None

  def setSpectrometerName(self, event):

    if self.spectrometer:
      default = 'eciSpectrometer%d' % self.spectrometer.serial
      self.spectrometer.name = self.spectrometerNameEntry.get() or default

  def setSpectrometerFreq(self, event):

    if self.spectrometer:
      value = self.spectrometerFreqEntry.get() or None
      self.spectrometer.setProtonFreq( value )

      if value is not None:
        value = '%d' % round(value)

      self.spectrometer.setNominalFreq( value )
      self.updateInstrumentsAfter()

  def setSpectrometerManufacturer(self, event):

    if self.spectrometer:
      manuName = self.spectrometerManufacturerPulldown.getObject()

      if manuName is None:
        msg = 'Other spectrometer manufacturer'
        manuName = askString('Input', msg, parent=self) or ''
        manuName.strip()

      if manuName != '':
        affStore = getAffiliationStore(self.project)

        manu = affStore.findFirstOrganisation(name=manuName)
        if manu is None:
          manu = affStore.newOrganisation(name=manuName)

        self.spectrometer.setManufacturer( manu )

      else:
        self.spectrometer.manufacturer = None

      #self.updateInstrumentsAfter()

  def setSpectrometerModel(self, event):

    if self.spectrometer:
      #self.spectrometer.setModel( self.spectrometerModelEntry.get() or None )

      modelName = self.spectrometerModelPulldown.getObject()

      if modelName is None:
        msg = 'Other spectrometer model'
        modelName = askString('Input', msg, parent=self) or ''
        modelName.strip()

      if modelName == '':
        showWarning('Failure','Spectrometer model cannot be empty.', parent=self)
        return

      self.spectrometer.model = modelName

  def setSpectrometerSerial(self, event):

    if self.spectrometer:
      self.spectrometer.serialNumber = self.spectrometerSerialEntry.get() or None

  def setSpectrometerDetails(self, event):

    if self.spectrometer:
      self.spectrometer.setDetails( self.spectrometerDetailsEntry.get() or None )

  def getSpectrometerName(self, spectrometer):

    text = ''
    if spectrometer:
      text = spectrometer.name

    self.spectrometerNameEntry.set(text)

  def getSpectrometerFreq(self, spectrometer):

    value = 0.0
    if spectrometer:
      value = spectrometer.protonFreq

      if not spectrometer.protonFreq:
        if spectrometer.nominalFreq:
          spectrometer.protonFreq = float(spectrometer.nominalFreq)
          value = spectrometer.nominalFreq

    self.spectrometerFreqEntry.set(value)  

  def getSpectrometerManufacturer(self, spectrometer):

    if spectrometer:
      index = 0
      names = SPEC_MANU_LIST[:]
      manus = SPEC_MANU_LIST[:]

      names.append('<Custom>')
      manus.append(None)

      if spectrometer.manufacturer:
        manuName = spectrometer.manufacturer.name

        if manuName not in manus:
          names.insert(0, manuName)
          manus.insert(0, manuName)

        index = manus.index(manuName)

      self.spectrometerManufacturerPulldown.setup(names, manus, index)

  def getSpectrometerModel(self, spectrometer):

    if spectrometer:
      modelName = spectrometer.model

      index = 0
      names = []
      models = []

      if spectrometer.manufacturer:
        for (manu, model) in SPEC_MODEL_LIST:
          if manu == spectrometer.manufacturer.name:
            names.append(model)
            models.append(model)

      names.append('<Other>')
      models.append(None)

      if modelName != '':

        if modelName not in models:
          names.insert(0, modelName)
          models.insert(0, modelName)

        index = models.index(modelName)

      self.spectrometerModelPulldown.setup(names, models, index)

  def getSpectrometerSerial(self, spectrometer):

    text = ''
    if spectrometer:
      text = spectrometer.serialNumber

    self.spectrometerSerialEntry.set(text)

  def getSpectrometerDetails(self, spectrometer):

    text = ''
    if spectrometer:
      text = spectrometer.details

    self.spectrometerDetailsEntry.set(text)


  # Probe functions

  def selectProbe(self, obj, row, col):

    self.probe = obj

  def addProbe(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    instrumentStore = getInstrumentStore(self.project)

    name = ''
    while not name:
      name = askString('Text Entry','Enter NMR Probe Name')

    probe = instrumentStore.newNmrProbe(name=name)

  def removeProbe(self):

    if self.probe:
      msg = 'Really remove specification of NMR Probe %s?' % self.probe.name

      if showOkCancel('query',msg):
        self.probe.delete()
        self.probe = None

  def setProbeName(self, event):

    if self.probe:
      default = 'eciProbe%d' % self.probe.serial
      self.probe.name = self.probeNameEntry.get() or default

  def setProbeType(self, index, name=None):

    if name is None:
      index = self.probeTypePulldown.getSelectedIndex()

    if self.probe:
      self.probe.setProbeType( NMR_PROBE_TYPES[index] )
      self.updateInstrumentsAfter()

  def setProbeManufacturer(self, event):

    if self.probe:
      value = self.probeManufacturerEntry.get() or None

      if value is not None:
        affStore = getAffiliationStore(self.project)

        manu = affStore.findFirstOrganisation(name=value)
        if manu is None:
          manu = affStore.newOrganisation(name=value)

        self.probe.setManufacturer( manu )

      #self.updateInstrumentsAfter()

  def setProbeModel(self, event):

    if self.probe:
      self.probe.setModel( self.probeModelEntry.get() or None )

  def setProbeSerial(self, event):

    if self.probe:
      self.probe.serialNumber = self.probeSerialEntry.get() or None

  def setProbeDiameter(self, event):

    if self.probe:
      self.probe.setDiameter( self.probeDiameterEntry.get() or None )
      self.updateInstrumentsAfter()

  def setProbeDetails(self, event):

    if self.probe:
      self.probe.setDetails( self.probeDetailsEntry.get() or None )

  def getProbeName(self,probe):

    text = ''
    if probe:
      text = probe.name

    self.probeNameEntry.set(text)

  def getProbeType(self,probe):

    index = -1
    if probe and probe.probeType:
      index = NMR_PROBE_TYPES.index(probe.probeType)

    self.probeTypePulldown.setIndex(index)

  def getProbeManufacturer(self, probe):

    text = ''
    if probe:
      manu = probe.manufacturer
      if manu:
        text = manu.name

    self.probeManufacturerEntry.set(text)

  def getProbeModel(self,probe):

    text = ''
    if probe:
      text = probe.model

    self.probeModelEntry.set(text)

  def getProbeSerial(self,probe):

    text = ''
    if probe:
      text = probe.serialNumber

    self.probeSerialEntry.set(text)

  def getProbeDiameter(self,probe):

    value = 0.0
    if probe:
      value = probe.diameter

    self.probeDiameterEntry.set(value)

  def getProbeDetails(self,probe):
  
    text = ''
    if probe:
      text = probe.details

    self.probeDetailsEntry.set(text)

  def updateInstrumentsAfter(self, *opt):

    #if self.waiting:
    #  return
    #else:
    #  self.waiting = True
    
    self.after_idle(self.updateInstruments)

  def updateInstruments(self):

    instrumentStore = getInstrumentStore(self.project)

    if not instrumentStore:
      return

    getInstruments = instrumentStore.findAllInstruments

    # Probes

    objectList = []
    textMatrix = []

    probes = [(p.serial, p) for p in getInstruments(className='NmrProbe')]
    probes.sort()
    for serial, probe in probes:
      manu = probe.manufacturer
      probeManu = ''
      if manu:
        probeManu = manu.name

      datum = [serial,
               probe.name,
               probe.probeType,
               probeManu,
               probe.model,
               probe.serialNumber,
               probe.diameter,
               probe.details]
 
      objectList.append(probe)
      textMatrix.append(datum)

    self.probeMatrix.update(objectList=objectList, textMatrix=textMatrix)

    # Spectrometers

    objectList = []
    textMatrix = []
    spectrometers = [(s.serial, s) for s in getInstruments(className='NmrSpectrometer')]
    spectrometers.sort()
    for serial, spectrometer in spectrometers:
      manu = spectrometer.manufacturer
      specManu = ''
      if manu:
        specManu = manu.name

      datum = [serial,
               spectrometer.name,
               spectrometer.nominalFreq,
               spectrometer.protonFreq,
               specManu,
               spectrometer.model,
               spectrometer.serialNumber,
               spectrometer.details]

      objectList.append(spectrometer)
      textMatrix.append(datum)

    self.spectrometerMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waiting = False


  # Sample Condition Set functions

  def selectSampleConditionSet(self, sampleConditionSet, row, col):

    self.sampleConditionSet = sampleConditionSet

    self.updateSampleConditionsAfter()

  def updateSampleConditionSetsAfter(self, *opt):

    #if self.waiting:
    #  return
    #else:
    #  self.waiting = True
    
    self.after_idle(self.updateSampleConditionSets)

  def updateSampleConditionSets(self, obj=None):

    textMatrix  = []
    objectList  = []

    if self.project and self.nmrProject:  

      for scs in self.nmrProject.sortedSampleConditionSets():

        if not scs.name:
          scs.name = 'CondSet' + str(scs.serial)

        datum = [scs.serial,
                 scs.name,
                 scs.details]

        objectList.append(scs)
        textMatrix.append(datum)

    self.nmrCondSetTable.update(textMatrix=textMatrix,
                                objectList=objectList)

    self.waiting = False

  def addSampleConditionSet(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    if not self.nmrProject:
      self.nmrProject = getNmrProject(self.project)

    msg = 'Name for sample condition set'
    cName = askString('Query', msg, parent=self) or ''
    cName = cName.strip()

    if not cName:
      showWarning('Failure', 'No name specified.', parent=self)
      return

    if self.nmrProject.findFirstSampleConditionSet(name=cName):
      showWarning('Failure', 'Type already in use.', parent=self)
      return

    self.sampleConditionSet = self.nmrProject.newSampleConditionSet(name=cName)

  def removeSampleConditionSet(self):

    if self.sampleConditionSet:
      self.sampleConditionSet.delete()
      self.sampleConditionSet = None

  def setSampleConditionSet(self, event):

    if self.sampleConditionSet:
      value = self.sampleConditionSetEntry.get() or None

      self.sampleConditionSet.name = value

  def setSampleConditionSetDetails(self, event):

    if self.sampleConditionSet:
      value = self.sampleConditionSetDetailsEntry.get() or None

      if value != '':
        self.sampleConditionSet.details = value

  def getSampleConditionSet(self, sampleConditionSet):

    value = None
    if sampleConditionSet:
      value = sampleConditionSet.name

    self.sampleConditionSetEntry.set(value)

  def getSampleConditionSetDetails(self, sampleConditionSet):

    value = None
    if sampleConditionSet:
      value = sampleConditionSet.details

    self.sampleConditionSetDetailsEntry.set(value)


  # Sample Condition functions

  def selectSampleCondition(self, sampleCondition, row, col):

    self.sampleCondition = sampleCondition

  def updateSampleConditionsAfter(self, *opt):

    #if self.waiting:
    #  return
    #else:
    #  self.waiting = True
    
    self.after_idle(self.updateSampleConditions)

  def updateSampleConditions(self, obj=None):

    if not self.sampleConditionSet:
      return

    textMatrix  = []
    objectList  = []

    if self.project and self.nmrProject:  

      for sc in self.sampleConditionSet.sortedSampleConditions():

        datum = [sc.condition,
                 sc.value,
                 sc.error,
                 sc.unit]

        objectList.append(sc)
        textMatrix.append(datum)

    self.nmrCondTable.update(textMatrix=textMatrix,
                             objectList=objectList)

    self.waiting = False

  def addStdSampleCondition(self):

    if not self.sampleConditionSet:
      showWarning('Failure','No condition set selected in top frame.', parent=self)
      return

    stdSampleConditions = ( ('temperature', 0.0, 'K'),
                            ('pH', 0.0, 'pH'),
                            ('pressure', 1.0, 'atm'),
                            ('ionic strength', 0.0, 'M') )

    for (cType, cVal, cUnit) in stdSampleConditions:

      if self.sampleConditionSet.findFirstSampleCondition(condition=cType):
        continue

      self.sampleCondition = self.sampleConditionSet.newSampleCondition(condition=cType,
                                                                        value=cVal,
                                                                        unit=cUnit)

  def addSampleCondition(self):

    if not self.sampleConditionSet:
      showWarning('Failure','No condition set selected in top frame.', parent=self)
      return

    if not self.nmrProject:
      self.nmrProject = getNmrProject(self.project)

    # TODO: do this as a list?

    msg = 'Type of sample condition'
    cType = askString('Query', msg, parent=self) or ''
    cType.strip()

    if cType == '':
      showWarning('Failure','Type of sample condition cannot be empty.', parent=self)
      return

    if self.sampleConditionSet.findFirstSampleCondition(condition=cType):
      showWarning('Failure','Type already in use.', parent=self)
      return

    self.sampleCondition = self.sampleConditionSet.newSampleCondition(condition=cType)

  def removeSampleCondition(self):

    if self.sampleCondition:
      self.sampleCondition.delete()
      self.sampleCondition = None

  def setSampleCondition(self, event):

    if self.sampleCondition:
      cType = self.sampleConditionPulldown.getObject()

      if cType is None:
        msg = 'Non standard sample condition type:'
        cType = askString('Input', msg, parent=self) or ''
        cType.strip()

      self.sampleCondition.condition = cType or None

  def setSampleConditionValue(self, event):

    if self.sampleCondition:
      value = self.sampleConditionValueEntry.get() or 0.0

      self.sampleCondition.value = value

  def setSampleConditionError(self, event):

    if self.sampleCondition:
      value = self.sampleConditionErrorEntry.get() or None

      self.sampleCondition.error = value

  def setSampleConditionUnit(self, event):

    if self.sampleCondition:
      unit = self.sampleConditionUnitEntry.get() or None
      self.sampleCondition.unit = unit

  def getSampleCondition(self, sampleCondition):

    if not hasattr(self, 'cTypes'):
      self.cTypes = SAMPLE_CONDITION_TYPES
      self.cTypes.sort()
      self.cTypes.append(None)

    cType = None
    if self.sampleCondition:
      cType = self.sampleCondition.condition

    if cType:
      if cType not in self.cTypes:
        self.cTypes.insert(0, cType)
        self.cTypes.remove(None)
        self.cTypes.sort()
        self.cTypes.append(None)

    index = self.cTypes.index(cType)

    names = self.cTypes[:-1] + ['<Other>',]

    self.sampleConditionPulldown.setup(names, self.cTypes, index)

  def getSampleConditionValue(self, sampleCondition):

    value = None
    if sampleCondition:
      value = sampleCondition.value

    self.sampleConditionValueEntry.set(value)

  def getSampleConditionError(self, sampleCondition):

    value = None
    if sampleCondition:
      value = sampleCondition.error

    self.sampleConditionErrorEntry.set(value)

  def getSampleConditionUnit(self, sampleCondition):

    text = ''
    if sampleCondition:
      text = sampleCondition.unit

    self.sampleConditionUnitEntry.set(text)


  # Structure generation functions

  def selectStrucGen(self, strucGen, row, col):

    self.strucGen = strucGen

  def addStrucGen(self):

    if not self.entry:
      showWarning('Failure','You need to add an "Entry" object in the "Main" tab.', parent=self)
      return

    sgName = askString('Query','Structure generation name:',parent=self) or ''
    sgName.strip()

    if sgName:
      strucGen = self.nmrProject.findFirstStructureGeneration(name=sgName)
      if not strucGen:
        strucGen = self.nmrProject.newStructureGeneration(name=sgName)

      #self.entry.addAuthor(person)

  def removeStrucGen(self):

    if not self.entry:
      return

    if self.strucGen:

      msg = 'Delete structure generation with name %s?'
      if showYesNo('Query', msg % self.strucGen.name, parent=self):
        self.strucGen.delete()
        self.updateStrucGens()

  def toggleStrucGen(self, strucGen):

    if strucGen and self.entry:
      if strucGen in self.entry.structureGenerations:
        strucGen.removeEntry(self.entry)
      else:
        strucGen.addEntry(self.entry)

  def updateStrucGenAfter(self, strucGen):

    self.after_idle(self.updateStrucGens)

  def updateStrucGens(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    nCols = len(self.strucGenTable.headingList)

    if self.entry and self.nmrProject:
      #for i, author in enumerate(self.entry.authors):
      for i, strucGen in enumerate(self.nmrProject.sortedStructureGenerations() ):

        strucGenSoft = None
        if strucGen.method:
          if strucGen.method.software:
            strucGenSoft = strucGen.method.software.name

        if strucGen and (strucGen in self.entry.structureGenerations):
          colors = [NICE_GREEN] * nCols
          submit = 'Yes'
        else:
          colors = [None] * nCols
          submit = 'No'

        strucEns = strucGen.structureEnsemble

        ensTitle = None

        if strucEns:
          ensTitle = '%s : %d' % (strucEns.molSystem.code, strucEns.ensembleId)

        constraintSet = strucGen.nmrConstraintStore

        consTitle = None

        if constraintSet:
          consTitle = '%s : %d' % (constraintSet.nmrProject.name, constraintSet.serial)

        datum = [i+1,
                 strucGen.name,
                 submit,
                 strucGen.generationType,
                 strucGen.details,
                 strucGenSoft,
                 ensTitle,
                 consTitle]

        objectList.append(strucGen)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    self.strucGenTable.update(textMatrix=textMatrix,
                              objectList=objectList,
                              colorMatrix=colorMatrix)

    self.updateStrucEns()
    self.updateConstraintSetPulldown()

  def getStrucGenName(self, strucGen):

    self.strucGenNameEntry.set(strucGen.name)

  def getStrucGenType(self, strucGen):

    if strucGen.generationType not in GENERATION_TYPES:
      strucGen.generationType = GENERATION_TYPES[1]

    index = GENERATION_TYPES.index(strucGen.generationType)

    self.strucGenTypePulldown.setup(GENERATION_TYPES, GENERATION_TYPES, index)

  def getStrucGenDetails(self, strucGen):

    self.strucGenDetailsEntry.set(strucGen.details)

  def getStrucGenSoftware(self, strucGen):

    if strucGen:
      methodStore = getMethodStore(self.project)

      index = 0
      names = []
      softs = []

      for soft in methodStore.sortedSoftware():
        names.append(soft.name)
        softs.append(soft)

      names.append('<None>')
      softs.append(None)

      if strucGen.method and strucGen.method.software:
        index = softs.index(strucGen.method.software)

      self.strucGenSoftwarePulldown.setup(names, softs, index)

  def setStrucGenName(self, event):

    if self.strucGen:
      name = self.strucGenNameEntry.get()

      self.strucGen.name = name

  def setStrucGenType(self, event):

    if self.strucGen:
      self.strucGen.generationType = self.strucGenTypePulldown.getObject()

  def setStrucGenDetails(self, event):

    if self.strucGen:

      details = self.strucGenDetailsEntry.get()

      if details != '':
        self.strucGen.details = details

  def setStrucGenSoftware(self, event):

    if self.strucGen:
      software = self.strucGenSoftwarePulldown.getObject()

      if software:
        method = software.findFirstMethod()

        if not method:
          showWarning('Failure','Please set a method for this software in the References tab.', parent=self)
          self.tabbedFrame.select(2)
          return

        else:
          #self.strucGen.__dict__['method'] = method
          #self.updateNmrListSelectTable()

          self.strucGen.setMethod(method)

      else:
        self.strucGen.setMethod(None)

  def setStrucEns(self, event):

    strucEns = self.strucEnsPulldown.getObject()

    if strucEns is None:
      self.strucGen.structureEnsemble = None

    else:
      self.strucGen.structureEnsemble = strucEns

  def setStrucGenConstSet(self, event):

    constraintSet = self.strucGenConstSetPulldown.getObject()

    if constraintSet is None:
      self.strucGen.nmrConstraintStore = None

    else:
      self.strucGen.nmrConstraintStore = constraintSet

  def getStrucEns(self, strucGen):

    index = 0
    names = []
    strucEnss = []

    if self.molSystem:
      strucEnss.extend(self.molSystem.sortedStructureEnsembles() )

    if strucEnss:
      strucEns = strucGen.structureEnsemble

      if strucEns not in strucEnss:
        strucEns = strucEnss[0]

      names = [ '%s : %d' % (se.molSystem.code, se.ensembleId) for se in strucEnss ]
      index = strucEnss.index(strucEns)

    names.append('<None>')
    strucEnss.append(None)

    self.strucEnsPulldown.setup(names, strucEnss, index)

  def getStrucGenConstSet(self, strucGen):

    index = 0
    names = []
    constraintSets = self.project.sortedNmrConstraintStores()

    if constraintSets:
      if strucGen.nmrConstraintStore and strucGen.nmrConstraintStore not in constraintSets:
        strucGen.nmrConstraintStore = constraintSets[0]

      names = [ '%s : %d'  % (cs.nmrProject.name, cs.serial) for cs in constraintSets ]

      if strucGen.nmrConstraintStore:
        index = constraintSets.index(strucGen.nmrConstraintStore)

    names.append('<None>')
    constraintSets.append(None)

    self.strucGenConstSetPulldown.setup(names, constraintSets, index)


  # Structure ensemble functions

  def selectStrucEns(self, ensemble, row, col):

    self.ensemble = ensemble

  def updateStrucEnsAfter(self, strucEns):

    self.after_idle(self.updateStrucEns)

  def updateStrucEns(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    nCols = len(self.strucEnsTable.headingList)

    if self.molSystem and self.entry:

      for i, ensemble in enumerate(self.molSystem.sortedStructureEnsembles() ):

        strucGen = ensemble.structureGeneration

        strucEnsCalc = None
        appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='calculated')
        if appData:
          strucEnsCalc = appData.value

        strucEnsCriteria = None
        appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='criteria')
        if appData:
          strucEnsCriteria = appData.value

        strucEnsRepr = None
        appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='representative')
        if appData:
          strucEnsRepr = appData.value

        strucEnsReprCri = None
        appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='repr_criteria')
        if appData:
          strucEnsReprCri = appData.value

        if strucGen and (strucGen in self.entry.structureGenerations):
          colors = [NICE_GREEN] * nCols
          submit = 'Yes'
        else:
          colors = [None] * nCols
          submit = 'No'

        chains = ensemble.sortedCoordChains()
        nRes = 0
        for chain in chains:
          nRes += len(chain.residues)

        datum = [ensemble.ensembleId,
                 ensemble.details,
                 len(ensemble.models),
                 ','.join([ch.code for ch in chains]),
                 nRes,
                 strucEnsCalc,
                 strucEnsCriteria,
                 strucEnsRepr,
                 strucEnsReprCri]

        objectList.append(ensemble)
        textMatrix.append(datum)
        colorMatrix.append(colors)

      self.checkEntryConsistency()

    self.strucEnsTable.update(textMatrix=textMatrix,
                              objectList=objectList,
                              colorMatrix=colorMatrix)

  def setStrucEnsDetails(self, event):

    if self.ensemble:

      details = self.strucEnsDetailsEntry.get()

      if details != '':
        self.ensemble.details = details

  def getStrucEnsDetails(self, ensemble):

    self.strucEnsDetailsEntry.set(ensemble.details)

  def setStrucEnsCalc(self, event):

    if self.ensemble:
      appData = self.ensemble.findFirstApplicationData(application='nmrStar',keyword='calculated')

      if appData:
        for appData in self.ensemble.findAllApplicationData(application='nmrStar',keyword='calculated'):
          self.ensemble.removeApplicationData(appData)

      text = self.strucEnsCalcEntry.get() or ''
      text.strip()

      if text:
        keywds = {'application': 'nmrStar',
                  'keyword':     'calculated',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.ensemble.addApplicationData(appData)

  def getStrucEnsCalc(self, ensemble):

    if ensemble:
      appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='calculated')

      if appData:
        self.strucEnsCalcEntry.set(appData.value)

  def setStrucEnsCriteria(self, event):

    if self.ensemble:
      appData = self.ensemble.findFirstApplicationData(application='nmrStar',keyword='criteria')

      if appData:
        for appData in self.ensemble.findAllApplicationData(application='nmrStar',keyword='criteria'):
          self.ensemble.removeApplicationData(appData)

      text = self.strucEnsCriteriaEntry.get() or ''
      text.strip()

      if text:
        keywds = {'application': 'nmrStar',
                  'keyword':     'criteria',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.ensemble.addApplicationData(appData)

  def getStrucEnsCriteria(self, ensemble):

    if ensemble:
      appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='criteria')

      if appData:
        self.strucEnsCriteriaEntry.set(appData.value)

  def setStrucEnsRepr(self, event):

    if self.ensemble:
      appData = self.ensemble.findFirstApplicationData(application='nmrStar',keyword='representative')

      if appData:
        for appData in self.ensemble.findAllApplicationData(application='nmrStar',keyword='representative'):
          self.ensemble.removeApplicationData(appData)

      text = self.strucEnsReprEntry.get() or ''
      text.strip()

      if text:
        keywds = {'application': 'nmrStar',
                  'keyword':     'representative',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.ensemble.addApplicationData(appData)

  def getStrucEnsRepr(self, ensemble):

    if ensemble:
      appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='representative')

      if appData:
        self.strucEnsReprEntry.set(appData.value)

  def setStrucEnsReprCri(self, event):

    if self.ensemble:
      appData = self.ensemble.findFirstApplicationData(application='nmrStar',keyword='repr_criteria')

      if appData:
        for appData in self.ensemble.findAllApplicationData(application='nmrStar',keyword='repr_criteria'):
          self.ensemble.removeApplicationData(appData)

      text = self.strucEnsReprCriEntry.get() or ''
      text.strip()

      if text:
        keywds = {'application': 'nmrStar',
                  'keyword':     'repr_criteria',
                  'value':       text}

        appData = Implementation.AppDataString(**keywds)
        self.ensemble.addApplicationData(appData)

  def getStrucEnsReprCri(self, ensemble):

    if ensemble:
      appData = ensemble.findFirstApplicationData(application='nmrStar',keyword='repr_criteria')

      if appData:
        self.strucEnsReprCriEntry.set(appData.value)


  # Constraint functions

  def updateConstraintSetPulldown(self):

    if not self.entry:
      return

    names = []
    index = 0
    constSet = None

    if self.project:
      constSets = self.project.sortedNmrConstraintStores()
    else:
      constSets = []

    if self.strucGen:
      constSet = self.strucGen.nmrConstraintStore

    elif self.entry.findFirstStructureGeneration():
      constSet = self.entry.findFirstStructureGeneration().nmrConstraintStore

    if constSets:
      names = [ '%s : %d' % (cs.nmrProject.name, cs.serial) for cs in constSets ]

      if constSet and constSet not in constSets:
        constSet = constSets[0]

      if constSet in constSets:
        index = constSets.index(constSet)  

    else:
      constSet = None

    if constSet is not self.constSet:
      self.changeConstraintSet(constSet)

    self.constraintSetPulldown.setup(names, constSets, index)
    self.updateConstraintLists()

  def changeConstraintSet(self, constSet):

    if constSet is not self.constSet:
      self.constSet = constSet
      self.updateConstraintLists()

  def selectConstraintList(self, constList, row, col):

    self.constList = constList

  def updateConstraintListsAfter(self, constList):

    self.after_idle(self.updateConstraintLists)

  def updateConstraintLists(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    constSet = None

    if self.entry:
      for strucGen in self.entry.structureGenerations:
        if self.constSet and self.constSet == strucGen.nmrConstraintStore:
          constSet = self.constSet
          break

    self.nmrProject = getNmrProject(self.project)

    if not self.constSet and self.nmrProject:
      self.constSet = getattr(self.project, 'currentNmrConstraintStore') or \
                      getattr(self.project, 'findFirstNmrConstraintStore')(nmrProject=self.nmrProject) or \
                      getattr(self.project, 'newNmrConstraintStore')(nmrProject=self.nmrProject)

    nCols = len(self.strucEnsTable.headingList)

    if self.entry and self.constSet:
      for constList in self.constSet.sortedConstraintLists():

        if constSet is not None:
          colors = [NICE_GREEN] * nCols
        else:
          colors = [None] * nCols

        constType = None
        if constList.className in CONSTRAINT_NAME_DATA:
          constType = CONSTRAINT_NAME_DATA[constList.className]

        datum = [constList.serial,
                 constType,
                 constList.name,
                 constList.details,
                 len(constList.constraints)]

        objectList.append(constList)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    self.constraintListsTable.update(textMatrix=textMatrix,
                                     objectList=objectList,
                                     colorMatrix=colorMatrix)

  def setConstListName(self, event):

    if self.constList:

      name = self.constListNameEntry.get()

      if name != '':
        self.constList.name = name

  def getConstListName(self, constList):

    self.constListNameEntry.set(constList.name)

  def setConstListDetails(self, event):

    if self.constList:

      details = self.constListDetailsEntry.get()

      if details != '':
        self.constList.details = details

  def getConstListDetails(self, constList):

    self.constListDetailsEntry.set(constList.details)


  # Citation functions

  def updateCitations(self, obj=None):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    nCols = len(self.citationTable.headingList)
    entry = self.entry

    if entry:
      citationStore = getCitationStore(self.project)

      for citation in citationStore.sortedCitations():
        #if citation.className != 'JournalCitation':
        #  continue  

        if not entry.primaryCitation:
          entry.primaryCitation = citation

        authors = []
        for author in citation.authors:
          if author.givenName:
            initials = ''.join([author.givenName[0]] + list(author.middleInitials))
          else:
            initials = ''.join(list(author.middleInitials))
          authors.append('%s %s' % (author.familyName, initials))   

        if len(authors) > 2:
          authors = authors[0] + ' et al.'
        elif len(authors) == 2:
          authors = ' & '.join(authors)
        elif authors:
          authors = authors[0]
        else:
          authors = None  

        editors = []
        for editor in citation.editors:
          if editor.givenName:
            initials = ''.join([editor.givenName[0]] + list(editor.middleInitials))
          else:
            initials = ''.join(list(editor.middleInitials))
          editors.append('%s %s' % (editor.familyName, initials))   

        if len(editors) > 2:
          editors = editors[0] + ' et al.'
        elif len(editors) == 2:
          editors = ' & '.join(editors)
        elif editors:
          editors = editors[0]
        else:
          editors = None  

        if citation is self.entry.primaryCitation:
          submit = 'primary'
          colors = [NICE_GREEN] * nCols
        elif citation in entry.otherCitations:
          submit = 'other'
          colors = ['#B0D0B0'] * nCols
        else:
          submit = 'no'
          colors = [None] * nCols

        citeType = 'Journal'

        journalAbb = None
        issue = None
        volume = None

        bookTitle = None
        publisher = None
        bookSeries = None
        city = None

        country = None
        institution = None

        confTitle = None

        if citation.className == 'JournalCitation':
          journalAbb = citation.journalAbbreviation
          issue = citation.issue
          volume = citation.volume
          citeType = 'Journal'
        elif citation.className == 'BookCitation':
          volume = citation.volume
          bookTitle = citation.bookTitle
          bookSeries = citation.bookSeries
          publisher = citation.publisher
          city = citation.publisherCity
          citeType = 'Book'
        elif citation.className == 'ThesisCitation':
          city = citation.city
          country = citation.country
          institution = citation.institution
          citeType = 'Thesis'
        elif citation.className == 'ConferenceCitation':
          city = citation.city
          country = citation.country
          confTitle = citation.conferenceTitle
          citeType = 'Conference'

        title = citation.title

        if title and len(title) > 60:
          title = title[:57]+'...'

        datum = [citation.serial,
                 citeType,
                 submit,
                 authors,
                 editors,
                 citation.year,
                 title,
                 journalAbb,
                 issue,
                 volume,
                 citation.firstPage,
                 citation.lastPage,
                 citation.status,
                 citation.pubMedId,
                 citation.doi,
                 ', '.join(citation.keywords),
                 citation.details,
                 bookTitle,
                 bookSeries,
                 publisher,
                 city,
                 country,
                 institution,
                 confTitle]

        objectList.append(citation)
        textMatrix.append(datum)
        colorMatrix.append(colors)

    self.citationTable.update(textMatrix=textMatrix,
                              objectList=objectList,
                              colorMatrix=colorMatrix)

  def doCiteTableEditMarkExtraRules(self, obj, row, col):

    if (col == 5 or col == 9 or col == 10 or col == 11 or col == 13 or col == 14) and (obj.status != 'published'):
      return False
    elif (col == 8) and (obj.status != 'published') and (obj.className == 'JournalCitation'):
      return False
    elif (col > 16) and (obj.className == 'JournalCitation'):
      return False
    elif (col == 7 or col == 8) and (obj.className != 'JournalCitation'):
      return False
    elif (col == 9) and (obj.className not in ('BookCitation', 'JournalCitation') ):
      return False
    elif (col > 20) and (obj.className == 'BookCitation'):
      return False
    elif (col == 17 or col == 18 or col == 19) and (obj.className != 'BookCitation'):
      return False
    elif (col == 22) and (obj.className != 'ThesisCitation'):
      return False
    elif (col == 23) and (obj.className != 'ConferenceCitation'):
      return False

    return True

  def getCiteDeposition(self, citation):

    status = citation.status
    if status not in CITATION_DEPOSIT:
      status = CITATION_DEPOSIT[-1]

    index = CITATION_DEPOSIT.index(status)

    self.citeDepositionPulldown.setup(CITATION_DEPOSIT, CITATION_DEPOSIT, index)

  def getCiteAuthors(self, citation):

    affStore = getAffiliationStore(self.project)

    selected = []
    names = []
    for person in affStore.sortedPersons():
      names.append(getPersonName(person))
      if person in citation.authors:
        selected.append(getPersonName(person))

    names.append('<None>')

    self.citationAuthorMulti.set(selected, names)

  def getCiteEditors(self, citation):

    affStore = getAffiliationStore(self.project)

    selected = []
    names = []
    for person in affStore.sortedPersons():
      names.append(getPersonName(person))
      if person in citation.editors:
        selected.append(getPersonName(person))

    names.append('<None>')

    self.citationEditorMulti.set(selected, names)

  def getCiteYear(self, citation):

    year = time.localtime()[0]

    cYear = citation.year or year

    years = range(min(cYear,year-37),max(cYear,year+3))
    years.reverse()

    index = years.index(cYear)
    names = [str(y) for y in years]

    self.citeYearPulldown.setup(names, years, index)

  def getCiteTitle(self, citation):

    text = citation.title or ''
    self.citeTitleEntry.setText(text)

  def getCiteJournal(self, citation):

    if citation.className != 'JournalCitation':
      return

    if not hasattr(self, 'journals'):
      self.journals = JOURNALS
      self.journals.sort()

    cJournal = citation.journalAbbreviation
    if cJournal:
      if cJournal not in self.journals:
        self.journals.append(cJournal)
        self.journals.sort()

    categories = [j[0] for j in self.journals] + [None]
    
    if cJournal:
      index = self.journals.index(cJournal)
    else:
      index = 0 # len(self.journals)

    names = self.journals + ['<Other>',]

    self.citeJournalPulldown.setup(names, self.journals+[True], index, categories=categories)

  def getCiteIssue(self, citation):

    if citation.className != 'JournalCitation':
      return

    text = citation.issue or ''
    self.citeIssueEntry.set(text)

  def getCiteVolume(self, citation):

    if citation.className not in ('JournalCitation', 'BookCitation'):
      return

    text = citation.volume or ''
    self.citeVolumeEntry.set(text)

  def getCiteFirstPage(self, citation):

    page = citation.firstPage
    self.citeFirstPageEntry.set(page)

  def getCiteLastPage(self, citation):

    page = citation.lastPage
    self.citeLastPageEntry.set(page)

  def getCiteStatus(self, citation):

    status = citation.status
    if status not in CITATION_STATUS:
      status = CITATION_STATUS[-1]

    index = CITATION_STATUS.index(status)

    self.citeStatusPulldown.setup(CITATION_STATUS, CITATION_STATUS, index)

  def getCitePubMed(self, citation):

    text = citation.pubMedId or ''
    self.citePubMedEntry.set(text)

  def getCiteDoi(self, citation):

    text = citation.doi or ''
    self.citeDoiEntry.set(text)

  def getCiteKeywords(self, citation):

    self.citeKeywordsMulti.setValues(citation.keywords)

  def getCiteDetails(self, citation):

    text = citation.details or ''
    self.citeDetailsEntry.set(text)

  def getBookTitle(self, book):

    if book.className != 'BookCitation':
      return

    text = book.bookTitle or ''
    self.bookTitleText.setText(text)

  def getBookPublisher(self, book):

    if book.className != 'BookCitation':
      return

    text = book.publisher or ''
    self.bookPublisherEntry.set(text)

  def getBookSeries(self, book):

    if book.className != 'BookCitation':
      return

    text = book.bookSeries or ''
    self.bookSeriesEntry.set(text)

  def getCiteCity(self, citation):

    text = ''

    if citation.className == 'BookCitation':
      text = citation.publisherCity or ''
    elif citation.className in ('ThesisCitation', 'ConferenceCitation'):
      text = citation.city or ''
    else:
      return

    self.citeCityEntry.set(text)

  def getCiteCountry(self, citation):

    if citation.className not in ('ThesisCitation', 'ConferenceCitation'):
      return

    text = citation.country or ''
    self.citeCountryEntry.set(text)

  def getThesisInst(self, thesis):

    if thesis.className != 'ThesisCitation':
      return

    text = thesis.institution or ''
    self.thesisInstitutionEntry.set(text)
    
  def getConfTitle(self, conf):

    if conf.className != 'ConferenceCitation':
      return

    text = conf.conferenceTitle or ''
    self.confTitleText.setText(text)

  def setCiteDeposition(self, event):

    choice = self.citeDepositionPulldown.getObject()
    citation = self.citation
    entry = self.entry
    primary = entry.primaryCitation

    if choice == CITATION_DEPOSIT[0]:

      if primary:
        if citation is not primary:
          if primary not in entry.otherCitations:
            entry.addOtherCitation(primary)

          entry.primaryCitation = citation

      else:
        if citation in entry.otherCitations:
          entry.removeOtherCitation(citation)

        entry.primaryCitation = citation

    elif choice == CITATION_DEPOSIT[1]:
      if citation is primary:
        entry.primaryCitation = None

      if citation not in entry.otherCitations:
        entry.addOtherCitation(citation)  

    else:
      if citation is primary:
        entry.primaryCitation = None

      if citation in entry.otherCitations:
        entry.removeOtherCitation(citation)

  def setCiteAuthors(self, event):

    affStore = getAffiliationStore(self.project)

    names = []
    people = []
    for person in affStore.persons:
      names.append(getPersonName(person))
      people.append(person)

    useNames = self.citationAuthorMulti.get()
    usePeople = []

    for name in useNames:
      if name == '<None>':
        continue
      index = names.index(name)
      person = people[index]

      if person not in usePeople:
        usePeople.append(person)

    self.citation.setAuthors(usePeople)
    self.citationTable.keyPressEscape()

  def setCiteEditors(self, event):

    affStore = getAffiliationStore(self.project)

    names = []
    people = []
    for person in affStore.persons:
      names.append(getPersonName(person))
      people.append(person)

    useNames = self.citationEditorMulti.get()
    usePeople = []

    for name in useNames:
      if name == '<None>':
        continue
      index = names.index(name)
      person = people[index]

      if person not in usePeople:
        usePeople.append(person)

    self.citation.setEditors(usePeople)
    self.citationTable.keyPressEscape()

  def setCiteYear(self, event):

    self.citation.year = self.citeYearPulldown.getObject()

  def setCiteTitle(self, event):

    text = self.citeTitleEntry.getText().strip() or ''
    text = ' '.join(text.split('\n')).strip()

    self.citation.title = text or None

  def setCiteJournal(self, event):

    if self.citation.className != 'JournalCitation':
      return

    journal = self.citeJournalPulldown.getObject()

    if journal is True:
      msg = 'Abbrevated journal name (e.g. "J. Mol. Biol.")'
      journal = askString('Input', msg, parent=self) or ''
      journal.strip()

    self.citation.journalAbbreviation = journal or None

  def setCiteIssue(self, event):

    if self.citation.className != 'JournalCitation':
      return

    self.citation.issue = self.citeIssueEntry.get() or None

  def setCiteVolume(self, event):

    if self.citation.className not in ('JournalCitation', 'BookCitation'):
      return

    self.citation.volume = self.citeVolumeEntry.get() or None

  def setCiteFirstPage(self, event):

    page = self.citeFirstPageEntry.get() or 0 
    page = str(page)

    self.citation.firstPage = page or None

  def setCiteLastPage(self, event):

    page = self.citeLastPageEntry.get() or 0
    page = str(page)

    self.citation.lastPage = page or None

  def setCiteStatus(self, event):

    self.citation.status = self.citeStatusPulldown.getObject()

  def setCitePubMed(self, event):

    self.citation.pubMedId = self.citePubMedEntry.get() or None

  def setCiteDoi(self, event):

    self.citation.doi = self.citeDoiEntry.get() or None

  def setCiteKeywords(self, event):

    keywords = [w.strip() for w in self.citeKeywordsMulti.get() if w.strip()]
    self.citationTable.keyPressEscape()
    self.citation.keywords = keywords 

  def setCiteDetails(self, event):

    self.citation.details = self.citeDetailsEntry.get() or None

  def setBookTitle(self, event):

    if self.citation.className != 'BookCitation':
      return

    text = self.bookTitleText.getText().strip() or ''
    text = ' '.join(text.split('\n')).strip()

    self.citation.bookTitle = text or None

  def setBookPublisher(self, event):

    if self.citation.className != 'BookCitation':
      return

    self.citation.publisher = self.bookPublisherEntry.get() or None

  def setBookSeries(self, event):

    if self.citation.className != 'BookCitation':
      return

    self.citation.bookSeries = self.bookSeriesEntry.get() or None

  def setCiteCity(self, event):

    if self.citation.className == 'BookCitation':
      self.citation.publisherCity = self.citeCityEntry.get() or None
    elif self.citation.className in ('ThesisCitation', 'ConferenceCitation'):
      self.citation.city = self.citeCityEntry.get() or None
    else:
      return

  def setCiteCountry(self, event):

    if self.citation.className not in ('ThesisCitation', 'ConferenceCitation'):
      return

    self.citation.country = self.citeCountryEntry.get() or None

  def setThesisInst(self, event):

    if self.citation.className != 'ThesisCitation':
      return

    self.citation.institution = self.thesisInstitutionEntry.get() or None

  def setConfTitle(self, event):

    if self.citation.className != 'ConferenceCitation':
      return

    text = self.confTitleText.getText().strip() or ''
    text = ' '.join(text.split('\n')).strip()

    self.citation.conferenceTitle = text or None

  def selectCitation(self, citation, row, col):

    self.citation = citation

  def addCitation(self):

    if self.entry:
      citationStore = getCitationStore(self.project)

      citationStore.newJournalCitation() #year=time.localtime()[0],
                                         #status=CITATION_STATUS[0])

  def addJournal(self):

    if self.entry:
      citationStore = getCitationStore(self.project)

      citationStore.newJournalCitation() #year=time.localtime()[0],
                                         #status=CITATION_STATUS[0])

  def addBook(self):

    if self.entry:
      citationStore = getCitationStore(self.project)

      citationStore.newBookCitation() #year=time.localtime()[0],
                                      #status=CITATION_STATUS[0])

  def addThesis(self):

    if self.entry:
      citationStore = getCitationStore(self.project)

      citationStore.newThesisCitation() #year=time.localtime()[0],
                                        #status=CITATION_STATUS[0])

  def addConference(self):

    if self.entry:
      citationStore = getCitationStore(self.project)

      citationStore.newConferenceCitation() #year=time.localtime()[0],
                                            #status=CITATION_STATUS[0])

  def removeCitation(self):

    if self.citation and showOkCancel('Confirm','Delete citation?',parent=self):
      self.citation.delete()
      self.citation = None


  # Software functions

  def updateSoftwareAfter(self, software):

    self.after_idle(self.updateSoftware)

  def updateSoftware(self, obj=None):

    textMatrix  = []
    objectList  = []

    entry = self.entry

    if entry:
      methodStore = getMethodStore(self.project)

      for software in methodStore.sortedSoftware():

        methodName = None

        method = software.findFirstMethod()

        if method:
          methodName = method.name

        datum = [software.name,
                 software.version,
                 ', '.join(software.tasks),
                 software.vendorName,
                 software.vendorAddress,
                 software.vendorWebAddress,
                 software.details,
                 methodName]

        objectList.append(software)
        textMatrix.append(datum)

      self.checkEntryConsistency()

    self.softwareTable.update(textMatrix=textMatrix,
                              objectList=objectList)

  def selectSoftware(self, software, row, col):

    self.software = software  

  def addSoftware(self):

    msg  = 'Name of software:'
    name = askString('Input', msg, parent=self) or ''
    name.strip()

    if not name:
      return

    msg = 'Software version (e.g. 2.1):'
    version = askString('Input', msg, parent=self) or ''
    version.strip()

    if not version:
      return

    methodStore = getMethodStore(self.project)
    software = methodStore.findFirstSoftware(name=name, version=version)
    if not software:
      software = methodStore.newSoftware(name=name, version=version)

  def removeSoftware(self):

    software = self.softwareTable.currentObjects

    if not software:
      return

    if len(software) == 1:
      msg = 'Delete software reference?'
    else:
      msg = 'Delete %s software references?' % len(software)

    if showOkCancel('Confirm', msg, parent=self):
      self.administerNotifiers(self.unregisterNotify)
      for softwareRef in software:
        softwareRef.delete()

      self.software = None  
      self.administerNotifiers(self.registerNotify)
      self.updateSoftware()

  def getSoftVersion(self, software):

    self.softVersionEntry.set(software.version)

  def getSoftTasks(self, software):

    self.softTasksMulti.set(software.tasks)

  def getSoftVendor(self, software):

    self.softVendorEntry.set(software.vendorName or '')

  def getSoftVendorAddr(self, software):

   self.softVendorAddrEntry.set(software.vendorAddress or '')

  def getSoftVendorWeb(self, software):

    self.softVendorWebEntry.set(software.vendorWebAddress or '')

  def getSoftDetails(self, software):

    self.softDetailsEntry.set(software.details or '')

  def getSoftMethod(self, software):

    if software.findFirstMethod():
      self.softMethodEntry.set(software.findFirstMethod().name)

  def setSoftVersion(self, event):

    text = self.softVersionEntry.get() or ''
    text.strip()

    self.software.version = text

  def setSoftTasks(self, event):

    tasks = self.softTasksMulti.get()
    tasks = [t.strip() for t in tasks]
    self.softwareTable.keyPressEscape()

    self.software.tasks = [t for t in tasks if t]

  def setSoftVendor(self, event):

    text = self.softVendorEntry.get() or ''
    text.strip()

    self.software.vendorName = text or None

  def setSoftVendorAddr(self, event):

    text = self.softVendorAddrEntry.get() or ''
    text.strip()

    self.software.vendorAddress = text or None

  def setSoftVendorWeb(self, event):

    text = self.softVendorWebEntry.get() or ''
    text.strip()

    self.software.vendorWebAddress = text or None

  def setSoftDetails(self, event):

    text = self.softDetailsEntry.get() or ''
    text.strip()

    self.software.details = text or None

  def setSoftMethod(self, event):

    text = self.softMethodEntry.get() or ''
    text.strip()

    if text != '':
      if not self.software.findFirstMethod(name=text):
        methodStore = getMethodStore(self.project)
        method = methodStore.newMethod(name=text)
        self.software.setMethods([method])

  def packageProject(self, filePrefix=None):

    import tarfile

    userPath = getUserDataPath(self.project)

    if not filePrefix:
      filePrefix = userPath

    tarFile = '%s.tgz' % filePrefix
    tarFp = tarfile.open(tarFile, 'w:gz')

    def visitDir(directory):
      files = os.listdir(directory)
      for relfile in files:
        fullfile = os.path.join(directory, relfile)
        if os.path.isdir(fullfile):
          visitDir(fullfile)
        elif relfile.endswith('.xml'):
          tarFp.add(fullfile)

    userDir = os.path.dirname(userPath)
    cwd = os.getcwd()
    os.chdir(userDir)
    try:
      userPath = os.path.basename(userPath)
      visitDir(userPath)
    finally:
      os.chdir(cwd)
      tarFp.close()

  def packageProject2(self):
    repository = self.project.findFirstRepository(name='userData')
    if not repository:
      return

    projectPath = repository.url.path
    parentDir, projectDir = os.path.split(projectPath)

    if not os.path.exists(projectDir):
      msg = 'Directory %s does not exist - save your CCPN project first.'
      showWarning('Warning', msg % projectDir)
      return

    writeable = os.access(parentDir, os.W_OK)
    if not writeable:
      msg = 'Directory %s does not have write access required to write temporary file.'
      showWarning('Warning', msg % parentDir)
      return

    tarName = '%s.tgz' % self.project.name
    fileName = os.path.join(parentDir, tarName)

    while os.path.exists(fileName):
      msg = 'Filename %s exists. Overwrite?' % fileName
      if showYesNo('Query', msg): 
        break

      tarName = '_' + tarName
      fileName = os.path.join(parentDir, tarName)

    cwd = os.getcwd()
    os.chdir(parentDir)

    tarFileObj = tarfile.open(fileName, 'w:gz')
    tarFileObj.add(projectDir)
    tarFileObj.close()

    os.chdir(cwd)

    showWarning('Info', 'File %s saved as tgz format' % tarName)

    return fileName

  def exportNmrStar31(self):

    if not self.project:
      showWarning('Failure','You need a CCPN project.', parent=self)
      return

    entryStore = self.project.currentNmrEntryStore or self.project.findFirstNmrEntryStore()

    if not entryStore:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    entry = entryStore.findFirstEntry()

    if not entry:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    fileTypes = [ FileType('STAR', ['*.str']),
                  FileType('All', ['*'])]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Export NMR-STAR 3.1 file', dismiss_text='Cancel',
                        selected_file_must_exist=False, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if os.path.exists(fileName):
      if not showYesNo('Save project', 'File already exists. Overwrite?', parent=self):
        return

    if fileName:
      nmrStarExport = NmrStarExport(entry, nmrStarVersion='3.1')

      nmrStarExport.createFile(fileName)
      nmrStarExport.writeFile()

  def exportPdb(self):

    if not self.project:
      showWarning('Failure','You need a CCPN project.', parent=self)
      return

    entryStore = self.project.currentNmrEntryStore or self.project.findFirstNmrEntryStore()

    if not entryStore:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    entry = entryStore.findFirstEntry()

    if not entry:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    nmrProject = self.project.currentNmrProject or self.project.findFirstNmrProject()

    if not nmrProject:
      showWarning('Failure','No NMR data in this project.', parent=self)
      return

    strucGen = entry.findFirstStructureGeneration()

    if not strucGen:
      showWarning('Failure','Please select a structure generation in the Structures tab.', parent=self)
      self.tabbedFrame.select(11)
      return

    strucEns = strucGen.structureEnsemble

    structures = []

    if strucEns:
      structures = strucEns.sortedModels()

    if not structures:
      showWarning('Failure','No structure models in this structure ensemble.', parent=self)
      self.tabbedFrame.select(11)
      return

    fileTypes = [ FileType('PDB', ['*.pdb']),
                  FileType('PDB Entry', ['*.ent']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Export PDB 3.20 file', dismiss_text='Cancel',
                        selected_file_must_exist=False, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if os.path.exists(fileName):
      if not showYesNo('Save project', 'File already exists. Overwrite?', parent=self):
        return

    if fileName:
      pdbObj = PdbFormat(self.project, self)

      pdbObj.writeCoordinates(fileName,
                              structures=structures,
                              resetMapping=False,
                              minimalPrompts=True,
                              verbose=True,
                              addPdbHeader=True,
                              version='3.20')


  def exportCnsDistance(self):

    if not self.project:
      showWarning('Failure','You need a CCPN project.', parent=self)
      return

    entryStore = self.project.currentNmrEntryStore or self.project.findFirstNmrEntryStore()

    if not entryStore:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    entry = entryStore.findFirstEntry()

    if not entry:
      showWarning('Failure','Please create a new CCPN Entry in the main tab.', parent=self)
      self.tabbedFrame.select(0)
      return

    nmrProject = self.project.currentNmrProject or self.project.findFirstNmrProject()

    if not nmrProject:
      showWarning('Failure','No NMR data in this project.', parent=self)
      return

    strucGen = entry.findFirstStructureGeneration()

    if not strucGen:
      showWarning('Failure','Please select a structure generation in the Structures tab.', parent=self)
      self.tabbedFrame.select(11)
      return

    nmrConstStore = strucGen.nmrConstraintStore
    distConstList = None
    distConsts = []

    if nmrConstStore:
      distConstList = nmrConstStore.findFirstConstraintList(className="DistanceConstraintList")

      if distConstList:
        distConsts = distConstList.sortedConstraints()

    if not distConsts:
      showWarning('Failure','No distance constraitns in this structure generation.', parent=self)
      self.tabbedFrame.select(11)
      return

    fileTypes = [ FileType('CNS', ['*.tbl']),
                  FileType('All', ['*']) ]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title='Export CNS file', dismiss_text='Cancel',
                        selected_file_must_exist=False, multiSelect=False,)

    fileName = fileSelectPopup.getFile()

    if os.path.exists(fileName):
      if not showYesNo('Save project', 'File already exists. Overwrite?', parent=self):
        return

    if fileName:
      cnsObj = CnsFormat(self.project, self)

      cnsObj.writeDistanceConstraints(fileName,
                                      constraintList=distConstList,
                                      resetMapping=False,
                                      minimalPrompts=True,
                                      verbose=True)


# General functions to be moved to a general spot later.

def getTopLevelStore(memopsRoot, className):

  if not memopsRoot:
    return

  store = getattr(memopsRoot, 'current'+className) or \
          getattr(memopsRoot, 'findFirst'+className)() or \
          getattr(memopsRoot, 'new'+className)(name='eciDefault')

  return store

def getTopLevelStore2(memopsRoot, className):

  if not memopsRoot:
    return

  store = getattr(memopsRoot, 'current'+className) or \
          getattr(memopsRoot, 'findFirst'+className)() or \
          getattr(memopsRoot, 'new'+className)(namingSystem='eciDefault')

  return store

def getMethodStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'MethodStore')

def getCitationStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'CitationStore')

def getAffiliationStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'AffiliationStore')

def getSampleStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'SampleStore')

def getRefSampleComponentStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'RefSampleComponentStore')

def getNmrProject(memopsRoot):

  return getTopLevelStore(memopsRoot, 'NmrProject')

def getEntryStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'NmrEntryStore') 

def getTaxonomy(memopsRoot):

  return getTopLevelStore(memopsRoot, 'Taxonomy') 

def getInstrumentStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'InstrumentStore') 

def getClassification(memopsRoot):

  return getTopLevelStore2(memopsRoot, 'Classification')

def getStructureGeneration(ensemble):

  strucGen = ensemble.structureGeneration

  if not strucGen:
    nmrProject = getNmrProject(ensemble.memopsRoot)
    strucGen = nmrProject.findFirstStructureGeneration()

    if not strucGen:
      strucGen = getNewStructureGeneration(ensemble, nmrProject)

    strucGen.structureEnsemble = ensemble

  return strucGen

def getNewStructureGeneration(ensemble, nmrProject):

  strucGen = nmrProject.newStructureGeneration(generationType=GENERATION_TYPES[1],
                                               name='eciDefault')

  strucGen.structureEnsemble = ensemble

  return strucGen

def setPersonInGroup(person, group):

  pinG = None

  if not person:
    return

  if group:
    pinG = person.currentPersonInGroup = person.newPersonInGroup(group=group)

  else:
    person.currentPersonInGroup = None

  return pinG

def getPersonInGroup(person):

  pinG = None

  if not person:
    return

  pinG = person.currentPersonInGroup

  # TODO: causes a bug when you remove a group
  #if not pinG:
  #  pinG = person.findFirstPersonInGroup()

  if pinG:
    person.currentPersonInGroup = pinG

  return pinG

def getExpSource(entryMolecule, orgName=UNDEFINED):

  if entryMolecule.experimentalSource:
    return entryMolecule.experimentalSource
  
  else:
    expSource = None
    entry = entryMolecule.entry
    
    """ 
    # wb104 8 Jul 2013: if the same source is used then editing the
    # Scientific Name of one, for example, means the same happens to all
    # the others, which is presumably normally not intended.

    # Use any existing sources first by default
    for entryMolecule2 in entry.entryMolecules:
      if entryMolecule2.experimentalSource:
        expSource = entryMolecule2.experimentalSource
        break
"""

    if not expSource:
      tax = getTaxonomy(entry.root)
      expSource = tax.findFirstNaturalSource(molecules=[]) or \
                  tax.newNaturalSource(scientificName=orgName, organismName=orgName)
    
    entryMolecule.experimentalSource = expSource
    
  return expSource

def getPersonName(person):

  data = (person.title,
          person.givenName,
          ' '.join(person.middleInitials),
          person.familyName,
          person.familyTitle)

  return ' '.join([d for d in data if d])

def getNmrLists(nmrProject, className, molSystem=None):

  if className == 'PeakList':
    peakLists = []

    for experiment in nmrProject.sortedExperiments():
      if molSystem and (molSystem not in experiment.molSystems):
        continue

      for spectrum in experiment.sortedDataSources():
        peakLists.extend(spectrum.sortedPeakLists())

    return peakLists

  elif className in DERIVED_LIST_CLASSES:
    derivedDataLists = []
    dLists = nmrProject.sortedDerivedDataLists()

    for dList in dLists:
      if dList.className == className:
        derivedDataLists.append(dList)

    return derivedDataLists

  else:
    measurementLists = []
    mLists = nmrProject.sortedMeasurementLists()

    for mList in mLists:
      if mList.className == className:
        measurementLists.append(mList)

    return measurementLists

def getNmrListName(nmrList):

  if nmrList.className == 'PeakList':
    spectrum = nmrList.dataSource
    experiment = spectrum.experiment

    return '%s : %s : %d' % (experiment.name, spectrum.name, nmrList.serial)

  else:
    return '%s : %s : %d' % (nmrList.className, nmrList.name, nmrList.serial)

#def analyseChemicalShifts(input):
#  return []

#def checkNmrEntryCompleteness(entry):
#  return [('Shift List', 'You must have at least one chemical shift list selected!','#FF4040'),]
