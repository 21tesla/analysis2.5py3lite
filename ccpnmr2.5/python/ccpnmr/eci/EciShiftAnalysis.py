from math import sqrt

PROTEIN_MOLTYPE = 'protein'
ROOT_TWO_PI = 2.506628274631

# Dicts to hel speed
REF_STORE_DICT = {}
CHEM_ATOM_REF_DICT = {}

REFDB_SD_MEAN =  {
  ('protein', 'Ala') : {
       "C":(177.828848,2.109680,0.452147,None),
       "CA":(53.352984,1.902676,0.210713,"HA"),
       "CB":(19.222672,1.746981,0.258369,"HB1"),
       "H":(8.198200,0.628909,0.033675,"N"),
       "HA":(4.276693,0.446260,0.146317,"CA"),
       "HB*":(1.364080,0.226681,0.212091,"CB"),
       "N":(122.708700,3.549744,0.159512,"H2"),
   },
   ('protein', 'Arg') : {
       "C":(176.476635,1.985958,0.496077,None),
       "CA":(56.956118,2.286571,0.236270,"HA"),
       "CB":(30.902098,1.776631,0.317743,"HB2"),
       "CD":(43.264272,0.631000,0.568196,"HD2"),
       "CG":(27.332456,1.047844,0.581473,"HG2"),
       "CZ":(160.254836,3.230362,0.982800,None),
       "H":(8.251157,0.626805,0.038322,"N"),
       "HA":(4.312683,0.477117,0.133977,"CA"),
       "HB2":(1.796660,0.261361,0.246530,"CB"),
       "HB3":(1.770165,0.267594,0.301147,"CB"),
       "HD2":(3.124174,0.221769,0.361798,"CD"),
       "HD3":(3.108990,0.237850,0.419433,"CD"),
       "HE":(7.358268,0.581758,0.718769,"NE"),
       "HG2":(1.565723,0.256262,0.342788,"CG"),
       "HG3":(1.543494,0.272259,0.396198,"CG"),
       "HH11":(6.806578,0.394752,0.968618,"NH1"),
       "HH12":(6.804218,0.413009,0.973446,"NH1"),
       "HH21":(6.725168,0.354246,0.971032,"NH2"),
       "HH22":(6.729098,0.436935,0.971636,"NH2"),
       "N":(120.208274,3.814534,0.208509,"H2"),
       "NE":(91.403172,13.311577,0.865419,"HE"),
       "NH1":(76.300233,12.184632,0.993060,"HH12"),
       "NH2":(75.668286,10.259347,0.993060,"HH21"),
   },
   ('protein', 'Asn') : {
       "C":(175.387335,1.722342,0.461435,None),
       "CA":(53.715068,1.780935,0.233736,"HA"),
       "CB":(38.882610,1.637450,0.291751,"HB2"),
       "CG":(176.826660,1.209325,0.911804,None),
       "H":(8.349260,0.632518,0.049966,"N"),
       "HA":(4.681485,0.371691,0.153253,"CA"),
       "HB2":(2.819383,0.300069,0.232059,"CB"),
       "HB3":(2.776859,0.314775,0.260899,"CB"),
       "HD21":(7.356185,0.463488,0.462441,"ND2"),
       "HD22":(7.141474,0.456559,0.466465,"ND2"),
       "N":(118.480741,4.131918,0.192153,"H2"),
       "ND2":(112.784090,2.263417,0.616030,"HD21"),
   },
   ('protein', 'Asp') : {
       "C":(176.552707,1.676911,0.428646,None),
       "CA":(54.879282,1.972805,0.195573,"HA"),
       "CB":(41.088030,1.475807,0.253906,"HB2"),
       "CG":(179.515246,1.549972,0.978906,None),
       "H":(8.315900,0.588151,0.036198,"N"),
       "HA":(4.603843,0.309474,0.133594,"CA"),
       "HB2":(2.726400,0.280057,0.220312,"CB"),
       "HB3":(2.683034,0.255266,0.255729,"CB"),
       "N":(120.230927,3.966225,0.146615,"H"),
   },
   ('protein', 'Cys') : { # CB from CSI, SD guess
       "C":(174.775094,1.963682,0.690521,None),
       "CA":(58.071390,3.292878,0.495920,"HA"),
       "CB":(28.6,1.8,0.549906,"HB2"),
       "H":(8.415828,0.683783,0.061519,"N"),
       "HA":(4.699157,0.549761,0.084118,"CA"),
       "HB2":(2.927206,0.473968,0.123038,"CB"),
       "HB3":(2.884519,0.495915,0.140615,"CB"),
       "HG":(2.539318,2.428652,0.998117,"SG"),
       "N":(119.339964,4.450665,0.438795,"H"),
   },
   ('protein', 'Cyss') : { # CB from CSI, SD guess
       "C":(174.775094,1.963682,0.690521,None),
       "CA":(58.071390,3.292878,0.495920,"HA"),
       "CB":(41.8,1.8,0.549906,"HB2"),
       "H":(8.415828,0.683783,0.061519,"N"),
       "HA":(4.699157,0.549761,0.084118,"CA"),
       "HB2":(2.927206,0.473968,0.123038,"CB"),
       "HB3":(2.884519,0.495915,0.140615,"CB"),
       "N":(119.339964,4.450665,0.438795,"H"),
   },
   ('protein', 'Gln') : {
       "C":(176.404030,1.920946,0.433155,None),
       "CA":(56.761798,2.082327,0.192157,"HA"),
       "CB":(29.358559,1.772749,0.260606,"HB3"),
       "CD":(179.827346,1.026666,0.926916,None),
       "CG":(33.876419,0.890980,0.516221,"HG2"),
       "H":(8.216737,0.612440,0.037790,"N"),
       "HA":(4.281435,0.448260,0.145811,"CA"),
       "HB2":(2.056673,0.235446,0.253476,"CB"),
       "HB3":(2.032193,0.243198,0.301604,"CB"),
       "HE21":(7.243197,0.430273,0.512656,"NE2"),
       "HE22":(7.014516,0.403005,0.512656,"NE2"),
       "HG2":(2.327917,0.242786,0.313369,"CG"),
       "HG3":(2.302919,0.262972,0.378610,"CG"),
       "N":(119.290354,3.686866,0.153298,"H"),
       "NE2":(111.888680,1.758582,0.612478,"HE21"),
   },
   ('protein', 'Glu') : {
       "C":(177.019843,1.918170,0.413230,None),
       "CA":(57.531344,2.064791,0.179611,"HA"),
       "CB":(30.197538,1.660817,0.241365,"HB2"),
       "CD":(182.899624,1.786294,0.982834,None),
       "CG":(36.202941,0.965699,0.497593,"HG3"),
       "H":(8.343540,0.613833,0.030563,"N"),
       "HA":(4.266130,0.434106,0.131673,"CA"),
       "HB2":(2.031746,0.198885,0.230061,"CB"),
       "HB3":(2.009291,0.202896,0.284069,"CB"),
       "HG2":(2.279046,0.202965,0.311911,"CG"),
       "HG3":(2.260935,0.206668,0.364455,"CG"),
       "N":(120.173056,3.592825,0.139837,"H2"),
   },
   ('protein', 'Gly') : {
       "C":(174.017069,1.841177,0.473382,None),
       "CA":(45.562676,1.121685,0.235463,"HA3"),
       "H":(8.330610,0.673335,0.044840,"N"),
       "HA2":(3.988405,0.369611,0.150491,"CA"),
       "HA3":(3.912364,0.373510,0.191237,"CA"),
       "N":(109.068021,3.809180,0.190418,"H"),
   },
   ('protein', 'His') : {
       "C":(175.302898,1.948641,0.514085,None),
       "CA":(56.675864,2.330521,0.279577,"HA"),
       "CB":(30.475373,2.009559,0.331690,"HB2"),
       "CD2":(119.952558,2.769141,0.787324,"HD2"),
       "CE1":(137.238155,2.346656,0.840141,"HE1"),
       "CG":(132.303333,3.787014,0.992253,None),
       "H":(8.262735,0.688922,0.100704,"N"),
       "HA":(4.627289,0.467128,0.186620,"CA"),
       "HB2":(3.110051,0.349972,0.261972,"CB"),
       "HB3":(3.057166,0.355392,0.291549,"CB"),
       "HD1":(8.265318,2.399673,0.943662,"ND1"),
       "HD2":(7.037220,0.429625,0.454930,"CD2"),
       "HE1":(8.013054,0.520085,0.517606,"CE1"),
       "HE2":(9.920824,2.588519,0.981690,"NE2"),
       "N":(119.109795,4.192319,0.235211,"H"),
       "ND1":(195.537424,37.008263,0.947183,"HD1"),
       "NE2":(177.698489,16.243778,0.952817,"HE2"),
   },
   ('protein', 'Ile') : {
       "C":(175.864340,1.828210,0.422507,None),
       "CA":(61.705518,2.600251,0.188136,"HA"),
       "CB":(38.854432,1.944834,0.252690,"HB"),
       "CD1":(13.584336,1.626539,0.482408,"HD11"),
       "CG1":(27.780248,1.684249,0.515848,"HG13"),
       "CG2":(17.636566,1.237755,0.469613,"HG22"),
       "H":(8.308229,0.703239,0.032858,"N"),
       "HA":(4.224109,0.580233,0.132306,"CA"),
       "HB":(1.784957,0.304183,0.218377,"CB"),
       "HD1*":(0.685842,0.282454,0.294272,"CD1"),
       "HG12":(1.270219,0.392634,0.329165,"CG1"),
       "HG13":(1.216866,0.399083,0.364059,"CG1"),
       "HG2*":(0.789078,0.258937,0.281477,"CG2"),
       "N":(121.173088,4.413042,0.150044,"H"),
   },
   ('protein', 'Leu') : {
       "C":(177.056361,1.919039,0.431646,None),
       "CA":(55.752106,2.056969,0.191568,"HA"),
       "CB":(42.540153,1.789209,0.255424,"HB2"),
       "CD1":(24.785797,1.489808,0.494267,"HD11"),
       "CD2":(24.227473,1.601505,0.521609,"HD23"),
       "CG":(26.875515,0.996005,0.543658,"HG"),
       "H":(8.238177,0.674706,0.033163,"N"),
       "HA":(4.346146,0.491254,0.154877,"CA"),
       "HB2":(1.616793,0.329304,0.261069,"CB"),
       "HB3":(1.534532,0.346944,0.302699,"CB"),
       "HD1*":(0.755018,0.268818,0.295114,"CD1"),
       "HD2*":(0.728231,0.285852,0.323337,"CD2"),
       "HG":(1.509274,0.327438,0.357559,"CG"),
       "N":(121.455301,4.025737,0.146587,"H2"),
   },
   ('protein', 'Lys') : {
       "C":(176.775235,1.906749,0.452691,None),
       "CA":(57.154802,2.135196,0.229246,"HA"),
       "CB":(32.992345,1.700470,0.286857,"HB3"),
       "CD":(29.036786,0.920116,0.588718,"HD2"),
       "CE":(41.966911,0.532535,0.611722,"HE2"),
       "CG":(25.002996,0.954247,0.555511,"HG3"),
       "H":(8.190187,0.626470,0.036807,"N"),
       "HA":(4.283917,0.449233,0.134827,"CA"),
       "HB2":(1.788608,0.231602,0.230646,"CB"),
       "HB3":(1.765164,0.240192,0.286457,"CB"),
       "HD2":(1.606115,0.230201,0.433887,"CD"),
       "HD3":(1.597734,0.235158,0.490898,"CD"),
       "HE2":(2.924529,0.170173,0.447289,"CE"),
       "HE3":(2.919496,0.175103,0.516303,"CE"),
       "HG2":(1.380507,0.240617,0.340068,"CG"),
       "HG3":(1.370381,0.245205,0.400480,"CG"),
       "HZ*":(7.483198,0.367244,0.962393,"NZ"),
       "N":(120.486754,3.843313,0.179236,"H"),
       "NZ":(48.671163,36.303503,0.998200,"HZ2"),
   },
   ('protein', 'Met') : {
       "C":(176.273626,2.049059,0.449133,None),
       "CA":(56.283714,2.173536,0.205727,"HA"),
       "CB":(33.259592,2.202185,0.269028,"HB2"),
       "CE":(17.113516,1.142950,0.691786,"HE1"),
       "CG":(32.098797,1.057559,0.567445,"HG3"),
       "H":(8.266266,0.624249,0.094951,"N"),
       "HA":(4.423979,0.489446,0.147702,"CA"),
       "HB2":(2.033462,0.323993,0.271289,"CB"),
       "HB3":(2.006533,0.336088,0.324039,"CB"),
       "HE*":(1.870402,0.366566,0.595328,"CE"),
       "HG2":(2.416216,0.370668,0.368500,"CG"),
       "HG3":(2.376700,0.421305,0.403919,"CG"),
       "N":(119.593417,3.693716,0.196684,"H2"),
   },
   ('protein', 'Phe') : {
       "C":(175.565521,1.948348,0.445665,None),
       "CA":(58.305101,2.551980,0.217338,"HA"),
       "CB":(40.181708,1.952764,0.279609,"HB2"),
       "CD1":(131.626637,1.261599,0.746439,"HD1"),
       "CD2":(131.620012,1.257999,0.828246,"HD2"),
       "CE1":(130.696973,1.484393,0.782662,"HE1"),
       "CE2":(130.731140,1.294712,0.853480,"HE2"),
       "CD*":(131.620012,1.257999,0.828246,"HD*"),
       "CE*":(130.696973,1.484393,0.782662,"HE*"),
       "CG":(137.410725,3.251523,0.994302,None),
       "CZ":(129.222986,1.622311,0.835165,"HZ"),
       "H":(8.382242,0.749659,0.058201,"N"),
       "HA":(4.648786,0.588854,0.160765,"CA"),
       "HB2":(2.998135,0.362851,0.245828,"CB"),
       "HB3":(2.952455,0.375443,0.268620,"CB"),
       "HD1":(7.058833,0.306262,0.365893,"CD1"),
       "HD2":(7.065550,0.304958,0.494506,"CD2"),
       "HE1":(7.081438,0.306502,0.435083,"CE1"),
       "HE2":(7.084839,0.303505,0.541718,"CE2"),
       "HD*":(7.065550,0.304958,0.494506,"CD*"),
       "HE*":(7.081438,0.306502,0.435083,"CE*"),
       "HZ":(6.996208,0.413472,0.582418,"CZ"),
       "N":(120.085399,4.232138,0.169312,"H3"),
   },
   ('protein', 'Pro') : {
       "C":(176.780325,1.499687,0.483755,None),
       "CA":(63.539193,1.404973,0.228159,"HA"),
       "CB":(31.977618,1.011604,0.292058,"HB2"),
       "CD":(50.448136,0.757834,0.560289,"HD2"),
       "CG":(27.325400,0.926700,0.575090,"HG3"),
       "HA":(4.408505,0.329217,0.180505,"CA"),
       "HB2":(2.077811,0.327969,0.257040,"CB"),
       "HB3":(2.029493,0.322236,0.275812,"CB"),
       "HD2":(3.665299,0.317100,0.331408,"CD"),
       "HD3":(3.643425,0.338100,0.362816,"CD"),
       "HG2":(1.944234,0.266628,0.361372,"CG"),
       "HG3":(1.920546,0.282568,0.404693,"CG"),
       "N":(131.085662,9.209093,0.979783,"H3"),
   },
   ('protein', 'Ser') : {
       "C":(174.683288,1.659248,0.464592,None),
       "CA":(58.850778,2.018005,0.217349,"HA"),
       "CB":(63.993613,1.384288,0.300172,"HB3"),
       "H":(8.290302,0.607958,0.058564,"N"),
       "HA":(4.513436,0.420120,0.147023,"CA"),
       "HB2":(3.884335,0.239598,0.248714,"CB"),
       "HB3":(3.861902,0.250070,0.313159,"CB"),
       "HG":(5.575657,1.195882,0.988238,"OG"),
       "N":(115.813074,3.723294,0.185739,"H3"),
   },
   ('protein', 'Thr') : {
       "C":(174.631688,1.707400,0.462333,None),
       "CA":(62.366001,2.580251,0.224368,"HA"),
       "CB":(69.886568,1.524338,0.292086,"HB"),
       "CG2":(21.649765,0.958760,0.519445,"HG23"),
       "H":(8.264453,0.646697,0.039162,"N"),
       "HA":(4.487641,0.495678,0.153658,"CA"),
       "HB":(4.167755,0.343464,0.253196,"CB"),
       "HG1":(5.132804,1.738798,0.975524,"OG1"),
       "HG2*":(1.144699,0.201156,0.278216,"CG2"),
       "N":(115.081474,4.991324,0.158009,"H"),
   },
   ('protein', 'Trp') : {
       "C":(176.218157,1.939642,0.518771,None),
       "CA":(57.833198,2.507920,0.324232,"HA"),
       "CB":(30.290294,1.875868,0.369738,"HB2"),
       "CD1":(126.600126,1.837343,0.715586,"HD1"),
       "CD2":(128.065000,2.873504,0.990899,None),
       "CE2":(137.984074,9.322232,0.980660,None),
       "CE3":(120.528840,1.591558,0.786121,"HE3"),
       "CG":(111.433182,0.972309,0.987486,None),
       "CH2":(123.902050,1.477779,0.773606,"HH2"),
       "CZ2":(114.407511,1.382949,0.754266,"HZ2"),
       "CZ3":(121.530652,1.444072,0.781570,"HZ3"),
       "H":(8.287413,0.814690,0.110353,"N"),
       "HA":(4.708245,0.556373,0.196815,"CA"),
       "HB2":(3.173023,0.343723,0.266212,"CB"),
       "HB3":(3.134429,0.352943,0.299204,"CB"),
       "HD1":(7.138713,0.340667,0.357224,"CD1"),
       "HE1":(10.116375,0.537535,0.356086,"NE1"),
       "HE3":(7.305717,0.382677,0.411832,"CE3"),
       "HH2":(6.964333,0.346475,0.420933,"CH2"),
       "HZ2":(7.286362,0.324283,0.379977,"CZ2"),
       "HZ3":(6.871129,0.368175,0.431172,"CZ3"),
       "N":(121.285471,4.389069,0.271900,"H2"),
       "NE1":(129.398739,2.021644,0.547213,"HE1"),
   },
   ('protein', 'Tyr') : {
       "C":(175.491657,1.900925,0.487770,None),
       "CA":(58.281804,2.475766,0.264346,"HA"),
       "CB":(39.616805,2.081812,0.337723,"HB3"),
       "CD1":(132.791984,1.567221,0.745532,"HD1"),
       "CD2":(132.624189,2.080215,0.835842,"HD2"),
       "CE1":(118.057560,1.375306,0.740828,"HE1"),
       "CE2":(118.004726,1.152713,0.833960,"HE2"),
       "CD*":(132.624189,2.080215,0.835842,"HD*"),
       "CE*":(118.057560,1.375306,0.740828,"HE*"),
       "CG":(129.985070,2.882228,0.988241,None),
       "CZ":(157.562708,1.371867,0.992004,None),
       "H":(8.325909,0.750383,0.062088,"N"),
       "HA":(4.645483,0.583003,0.161336,"CA"),
       "HB2":(2.903203,0.371810,0.248824,"CB"),
       "HB3":(2.849544,0.377313,0.265757,"CB"),
       "HD1":(6.937622,0.279576,0.331609,"CD1"),
       "HD2":(6.935375,0.282040,0.447789,"CD2"),
       "HE1":(6.704741,0.220055,0.355597,"CE1"),
       "HE2":(6.705084,0.218493,0.468956,"CE2"),
       "HD*":(6.935375,0.282040,0.447789,"CD*"),
       "HE*":(6.704741,0.220055,0.355597,"CE*"),
       "HH":(9.149362,1.563879,0.984478,"OH"),
       "N":(120.200602,4.345637,0.224365,"H"),
   },
   ('protein', 'Val') : {
       "C":(175.732054,1.819574,0.437803,None),
       "CA":(62.631900,2.794401,0.199169,"HA"),
       "CB":(32.935144,1.675944,0.260559,"HB"),
       "CG1":(21.586062,1.243471,0.477960,"HG12"),
       "CG2":(21.437030,1.431185,0.509808,"HG21"),
       "H":(8.307943,0.704303,0.039465,"N"),
       "HA":(4.204042,0.596150,0.143780,"CA"),
       "HB":(1.992687,0.286555,0.228479,"CB"),
       "HG1*":(0.832713,0.242229,0.252943,"CG1"),
       "HG2*":(0.810046,0.266148,0.275790,"CG2"),
       "N":(120.724161,4.673536,0.144242,"H"),
   }}

def analyseChemicalShifts(shiftList):
  
  project = shiftList.root
  data = []
    
  for shift in shiftList.measurements:
    resonance = shift.resonance
    resName = '%s:%s%8.8s%s' % getResonanceAtomTuple(resonance)
  
    resonance  = shift.resonance
    deltaMax   = None
    nContribs  = 0
    bmrbMean   = None
    randomCoil = None
    dev        = None
    
    resonanceSet = resonance.resonanceSet
    if resonanceSet:
      atomSet = resonanceSet.findFirstAtomSet()
      residue = atomSet.findFirstAtom().residue
      ccpCode = residue.ccpCode
      molType = residue.molResidue.molType
      atomName = atomSet.name
      chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType=molType)
      if chemAtomNmrRef is None:
        atomName = atomSet.findFirstAtom().name
        chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType=molType)
      
      if chemAtomNmrRef:
        dev = lookupAtomDeviation(project, ccpCode, atomName, shift.value, molType)
        bmrbMean   = chemAtomNmrRef.meanValue
        randomCoil = chemAtomNmrRef.randomCoilValue
  
    datum = [resonance.serial,resonance.isotopeCode,resName,bmrbMean,randomCoil,
             dev,shift.value,shift.error]
    
    data.append([shift, datum])
    
  return data  


def getResonanceAtomTuple(resonance):
  """Descrn: Give a tupe of string identifiers for a resonance
             indicating, chain, residue and atomic assignment
     Inputs: Nmr.Resonance
     Output: Tuple of Words (MolSystem.Chain.code,
             MolSystem.Residue identifier, Nmr.Resonance identifier)
  """
  
  molSystemCode = ''
  chainCode = ''
  res   = ''

  resonanceSet = resonance.resonanceSet
  spinSystem = resonance.resonanceGroup

  name = '[%d]' % resonance.serial
  if resonanceSet:
    atomSets     = tuple(resonanceSet.atomSets)
    
    atomNames = []
    for atomSet in atomSets:
      for atom in atomSet.atoms:
        atomNames.append(atom.name)

    tryName = atomNames[0]
    end = ''
    for atomName in atomNames[1:]:
      size = len(atomName)
      i = 0
      n = min(size, len(tryName))
      while i < n and atomName[i] == tryName[i]:
        i +=1
 
      tryName = tryName[0:i]
      end = '*'
 
    name = tryName[:1] + tryName[1:].lower() + end
    
    if len(atomSets) > 1:
      name = name[:-1] + getAmbigProchiralLabel(resonance)
      if len(atomSets[0].atoms) > 1:
        name = name + '*'
    
  elif resonance.assignNames:
    assignNames = tuple(resonance.assignNames)
    if assignNames:
      N = len(assignNames)
      if N > 1:
        for i in range(N):
          name = assignNames[i] + name
          if i < (N-1):
            name = '|' + name
      else:
        name = resonance.assignNames[0] + name
  
  elif resonance.name and (resonance.name != 'r%d' % resonance.serial):
    name = '[%d:%s]' % (resonance.serial,resonance.name)
  

  if resonanceSet:
    residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    chain = residue.chain
    chainCode = chain.code
    molSystemCode = chain.molSystem.code
    res = str(residue.seqCode) + residue.ccpCode
    
  elif spinSystem:
    residue = spinSystem.residue
    ccpCode = spinSystem.ccpCode
    residueProbs = [rp for rp in spinSystem.residueProbs if rp.weight > 0.0]
 
    if residue:
      sysChain = residue.chain
      chainCode = sysChain.code
      molSystemCode = sysChain.molSystem.code
      res = str(residue.seqCode) + residue.ccpCode

    elif residueProbs:
      resTexts = []
      resSeqs = []
      resCodes = set()
      
      for residueProb in residueProbs:
        residue = residueProb.possibility
        seq = residue.seqCode
        code = residue.ccpCode
        
        resTexts.append('%d?%s' % (seq,code))
        resSeqs.append('%d?' % seq)
        resCodes.add(code)
      
      if len(resCodes) == 1:
        res = '/'.join(resSeqs) + resCodes.pop()
      else:
        res = '/'.join(resTexts)
      
    elif ccpCode:
      res  = '{%d}%s' % (spinSystem.serial,ccpCode)
 
    else:
      res  = '{%d}' % (spinSystem.serial)
   
  return (molSystemCode, chainCode, res, name)

  
def getAtomGaussianDeviation(ccpCode, atomName, shiftValue, molType=PROTEIN_MOLTYPE):

  shiftRefs = REFDB_SD_MEAN.get((molType, ccpCode))
    
  if not shiftRefs:
    return 

  stats = shiftRefs.get(atomName)
  if not stats:
    return

  mean, sd, pMissing, bound = stats
  d = shiftValue-mean
  e = d/sd   
  #p = exp(-0.5*e*e)/(sd*ROOT_TWO_PI)
  
  return abs(e)  

def lookupAtomDeviation(project, ccpCode, atomName, shiftValue,
                        molType=PROTEIN_MOLTYPE):

  if molType == PROTEIN_MOLTYPE:
    deviation = getAtomGaussianDeviation(ccpCode, atomName, shiftValue)
  
  if deviation:
    return deviation
    
  chemAtomNmrRef = getChemAtomNmrRef(project, atomName, ccpCode, molType)
 
  if not chemAtomNmrRef:
    return
 
  distribution  = chemAtomNmrRef.distribution
  refPoint      = chemAtomNmrRef.refPoint
  refValue      = chemAtomNmrRef.refValue
  valuePerPoint = chemAtomNmrRef.valuePerPoint
  
  mean = 0.0
  mean2 = 0.0
  n = 0.0
  for i, prob in enumerate(distribution):
    shift = refValue+(i-refPoint)*valuePerPoint
    value = prob*shift
    mean += value
    mean2 += value*shift
    n += prob
    
  if n <= 0.0:
    return
  
  mean /= n
  mean2 /= n
  sigma2= abs(mean2 - (mean * mean))
  sigma = sqrt(sigma2)
  
  return abs((shiftValue-mean)/sigma)

def getChemAtomNmrRef(project, atomName, ccpCode, molType=PROTEIN_MOLTYPE,
                      sourceName='RefDB'):

  atomKey = molType + ccpCode + atomName
  chemAtomNmrRef  = CHEM_ATOM_REF_DICT.get(atomKey)

  if not chemAtomNmrRef:
    key = molType + ccpCode
    getRefStore = project.findFirstNmrReferenceStore
    nmrRefStore = REF_STORE_DICT.get(key, getRefStore(molType=molType,ccpCode=ccpCode)) 
    
    if nmrRefStore:
      REF_STORE_DICT[key] = nmrRefStore
      chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName=sourceName)
      if chemCompNmrRef:
        chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',descriptor='any')
        
        if chemCompVarNmrRef:
          for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
            if atomName == chemAtomNmrRef1.name:
              chemAtomNmrRef = chemAtomNmrRef1
              CHEM_ATOM_REF_DICT[atomKey] = chemAtomNmrRef
              break
        else:
          return
      else:
        return
    else:
      return
 
  if not chemAtomNmrRef:
    atomName2 = atomName[:-1]
    for chemAtomNmrRef1 in chemCompVarNmrRef.chemAtomNmrRefs:
      if atomName2 == chemAtomNmrRef1.name:
        chemAtomNmrRef = chemAtomNmrRef1
        CHEM_ATOM_REF_DICT[atomKey] = chemAtomNmrRef
        break

  return chemAtomNmrRef


def getAmbigProchiralLabel(resonance):
  """Descrn: Deterimine if an ambigous prochiral resonance (non-stereospecifically assigned)
             Has an "a" label or a "b" label. "a" is reserved for the upfield proton and any
             other nulceus bound to it.
     Inputs: Nmr.Resonance
     Output: Character
  """

  letter = ''
  resonanceSet = resonance.resonanceSet
  
  if resonanceSet:
    if resonance.isotopeCode == '1H':
      data = []
      for resonance2 in resonanceSet.sortedResonances():
        if resonance2.shifts:
          data.append( ('%f%d' % (resonance2.findFirstShift().value,resonance2.serial),resonance2) )
        else:
          data.append( (resonance2.serial,resonance2) )
 
      data.sort()
      resonances = [x[1] for x in data]
      i = resonances.index(resonance)
      letter = chr(ord('a')+i)

    else:
      resonance2 = resonance.findFirstCovalentlyBound(isotopeCode='1H')
 
      if resonance2 and resonance2.resonanceSet and \
           (len(resonance2.resonanceSet.atomSets) > 1):
        letter = getAmbigProchiralLabel(resonance2)
        resonance2.onebond = resonance
      
      elif (len(resonanceSet.resonances) > 1) and (len(resonanceSet.atomSets) > 1):
        for resonance2 in resonanceSet.resonances:
          if resonance2 is not resonance:
            resonance3 = resonance2.findFirstCovalentlyBound()
            if resonance3 and resonance3.resonanceSet and \
                 (len(resonance3.resonanceSet.atomSets) > 1):
              letter = 'b'
            break
             
      if not letter:
        data = []
        for resonance2 in resonanceSet.resonances:
          if resonance2.shifts:
            data.append( (resonance2.findFirstShift().value,resonance2) )
          else:
            data.append( (resonance2.serial,resonance2) )
 
        data.sort()
        resonances = [x[1] for x in data]
        i = resonances.index(resonance)
        letter = chr(ord('a')+i)       

  return letter
