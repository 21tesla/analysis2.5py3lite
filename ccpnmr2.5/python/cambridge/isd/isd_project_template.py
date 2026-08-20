##
## The Inferential Structure Determination (ISD) software library
##
## Authors: Wolfgang Rieping and Michael Habeck
##
##          Copyright (C) Michael Habeck and Wolfgang Rieping
##
##          All rights reserved.
##
## NO WARRANTY. This library is provided 'as is' without warranty of any
## kind, expressed or implied, including, but not limited to the implied
## warranties of merchantability and fitness for a particular purpose or
## a warranty of non-infringement.
##
## Distribution of substantively modified versions of this module is
## prohibited without the explicit permission of the copyright holders.
##

###############################################################################
##                                                                           ##
## Please do not modify the following lines                                  ##
##                                                                           ##
##############################################################################

import sys, os, commands

modules = os.path.join(os.environ['ISD_ROOT'], 'src', 'py')

if not modules in sys.path:
    sys.path.insert(0, modules)

from Isd.setup import *

sim_data = []

sim = SimulationSpecs()
sim.data_sets = sim_data
sim.isd_version = 1.2

###############################################################################
##                                                                           ##
## ISD project file (version 1.1)                                            ##
##                                                                           ##
###############################################################################
##                                                                           ##
## To modify the setup of an ISD calculation, please edit the settings below.##
##                                                                           ##
## Notes:                                                                    ##
##                                                                           ##
## - The environment variable ISD_ROOT needs to be set to the                ##
##   name of the directory in which ISD has been installed.                  ##
##                                                                           ##
###############################################################################
##                                                                           ##
## A quick guide to Python                                                   ##
##                                                                           ##
## An ISD project file is a Python script that is executed in order to start ##
## a calculation.                                                            ##
##                                                                           ##
## - Python expressions are case sensitive.                                  ##
## - Indentation matters, do not use whitespaces at the beginnin of a line.  ##
## - Comments start with a hash symbol (#).                                  ##
## - Strings are defined with a single (') or double (") quotation mark.     ##
## - The decimal point in floating point numbers is a dot, not a comma       ##
## - The value of a boolean variable can be True or False (case sensitive!). ##
## - Python dictionaries map keys to values.                                 ##
##   Example: d = {'name': 'mr smith', 'age': 123}                           ##
##                                                                           ##
###############################################################################
##                                                                           ##
## General settings                                                          ##
##                                                                           ##
###############################################################################

## name:             (string) The name of the simulation
## working_path:     (string) The directory in which ISD saves the results.
## temp_path:        (string) A directory accessible from all machines.
## shared_temp_path: (boolean) Set to True if temp_path is accessible from all
##                   machines (s. host_list below).
## cns_executable:   (string) The executable of the program CNS.
##                   This is required only if the sequence is specified in the
##                   form of a SEQ file.
## fileroot:         (string) Output files start with the fileroot,
##                   such as isd_13.pdb.
## temperature:      (float) The temperature of the system in Kelvin;
##                   default: 300.
## naming_system:    (const) Naming system used for sequence PDB file and data
##                   files; either IUPAC or CNS.

sim.name               = 'ISD simulation'
sim.working_path       = './'
sim.temp_path          = '/tmp'
sim.shared_temp_path   = True
sim.cns_executable     = commands.getoutput('which cns_solve')
sim.fileroot           = 'isd_sim'
sim.temperature        = 300.
sim.naming_system      = IUPAC

###############################################################################
##                                                                           ##
## Generation of PDB files                                                   ##
##                                                                           ##
###############################################################################

## ensemble:          (const) Set to REPRESENTATIVE to write representative
##                    ensemble members to disk. If set to MOST_PROBABLE, only
##                    the most probable structures is saved.
## n:                 (integer) Number of ensemble members that are stored as
##                    PDB file
## residue_list:      (list) List of residue numbers that shall be written
##                    to the PDB files. This feature enables one to exclude
##                    certain residues from the (WhatIf etc) quality analysis.
##                    Default: None (PDB files contain all residues).
##                    Currently, the number of the first residue is 1.

sim.pdb_files.ensemble     = REPRESENTATIVE
sim.pdb_files.n            = 100
sim.pdb_files.residue_list = None

###############################################################################
##                                                                           ##
## CCPN data model                                                           ##
##                                                                           ##
###############################################################################

## project_filename:      (string) Full path of the CCPN project that shall
##                        be used to import or export data.
## nmr_project_name:      (string) Name of the "NmrProject" within a CCPN
##                        project under which all exported objects (such as
##                        structures) are stored.
## molecular_system_name: (string) Name of the "MolSystem" that describes
##                        the structures to be exported.
## project.key:           (string) Unique key that is used to retrieve an ISD
##                        project in a CCPN project. ISD projects are stored
##                        under the NmrProject specified in "nmr_project_name"
## project.enabled:       (boolean) Set to True to export the simulation
##                        settings to CCPN (when using the command line option
##                        --ccpn-export). Default: False
## ensemble.enabled:      (boolean) Set to True to export the conformational
##                        ensemble to a CCPN project. The number of structures
##                        to be exported can be specified in section pdb_files.
##                        Molecular structures are stored under NmrProject
##                        with name "nmr_project_name".

sim.ccpn.project_filename               = ''
sim.ccpn.nmr_project_name               = ''
sim.ccpn.export.molecular_system_name   = ''
sim.ccpn.export.project.key             = ''
sim.ccpn.export.project.enabled         = True
sim.ccpn.export.ensemble.enabled        = True

###############################################################################
##                                                                           ##
## Molecular system                                                          ##
##                                                                           ##
###############################################################################

## filename:             (string) Either a PDB file or a text file that
##                       contains the sequence in three-letter codes.
## format:               (const) Format of sequence file, PDB, SEQ or XML.
## key:                  (string) Identifier used to retrieve a molecular
##                       sequence from a CCPN project.
## first_residue_number: (integer) Sequence number of the first resdiue;
##                       only if sequence is specified in SEQ format.Default:1.
## initial_conformation: (const) Defines the initial conformation used for the
##                       calculation: EXTENDED or RANDOM. Set to AS_IS to use
##                       the conformation contained in the PDB file.
## exclude_hydrogens:    (boolean) Set to True to exclude hydrogens from the
##                       calculation of the nonbonded energy.

sim.molecule.filename             = ''
sim.molecule.format               = CCPN
sim.molecule.key                  = ''
sim.molecule.first_residue_number = 1
sim.molecule.initial_conformation = EXTENDED
sim.molecule.exclude_hydrogens    = False

###############################################################################
##                                                                           ##
## Replica-exchange Monte Carlo                                              ##
##                                                                           ##
###############################################################################

## n_replicas:    (integer) Number of replicas; default: 50
## n_samples:     (integer) Number of replica exchange Monte Carlo samples
## hmc.steps:     (integer) Number of Hybrid Monte Carlo (HMC) step per
##                replica exchange step; default: 10
## hmc.md.steps:  (integer) Number of MD steps per HMC step; default: 250
## hmc.stepsize:  (float) Time step of the HMC leapfrog integration
## hmc.adjust_stepsize: (boolean/integer) Flag to switch the automatic
##                adjustment of the HMC step sizes on (True) or off (False);
##                if set to an integer, the step size adjustment will be
##                switched on until the number of samples exceeds the number
##                specified; the remaining samples will be generated with fixed
##                step sizes, however adjusted to each heatbath.
## communication: (string) Communication method. Either 'shared' or 'pyro'.
## save_interval: (integer) Interval in which states, that are sampled from the
##                posterior distribution, are saved.
## full_save:     (boolean) If set to True, all ensembles are written do disk
##                when a simulation is saved. Default: False
## python:        (string) Name of Python binary used on remote machines;
##                default: 'python'
## niceness:      (integer) Niceness of remote jobs; default: 0
## background:    (boolean) Set this variable to True if you wish to run ISD
##                without a terminal window (i.e. without interactive Python
##                shell). Default: False
## host_list:     (list of strings) Names or IP addresses of machines used for
##                the calculation. Names need to be specified as a
##                comma-separated list of strings, e.g. 'node1', '192.168.0.1'
## name_server:   (string) Name or IP address of the machine that runs the
##                PyRo name server.
## temp_paths     (dictionary) If sim.temp_path is not accessible from all
##                machines, you can specify the temporary directory
##                individually for each machine. If set to {}, it is assumed
##                that sim.temp_path exists locally on each machine.
##
## override_parameters: (boolean) When a calculation is resumed, set to True to
##                      override values of auxiliary parameters stored on disk
##                      with the ones specified in the project file. The
##                      variable is set to False by default, i.e. if restarted,
##                      a calculation continues where it stopped before.

sim.replica.n_replicas          = 50
sim.replica.n_samples           = 10000

sim.replica.hmc.steps           = 10
sim.replica.hmc.md.steps        = 250
sim.replica.hmc.stepsize        = 1e-3
sim.replica.hmc.adjust_stepsize = 1000

sim.replica.communication       = 'shared'
sim.replica.override_parameters = False
sim.replica.save_interval       = 50
sim.replica.full_save           = False
sim.replica.python              = commands.getoutput('which python')
sim.replica.niceness            = 10
sim.replica.background          = False

sim.replica.name_server         = 'localhost'
sim.replica.host_list           = 'localhost',
sim.replica.temp_paths          = {}

## Likelihood weight; acts similar to an inverse temperature.
##
## initial: (float) Initial likelihood weight; value must lie in interval
##          [0., 1]; dafault: 1.0
## final:   (float) Final likelihood weight; default: 0.05
## slope:   (float) Slope of power law used to calculate intermediate weights.
## first:   (integer) Index of the first replica for which the weight begins
##          to change. Weights of replicas with indices smaller than this index
##          are set to sim.weight_schedule.initial
## last:    (integer) Index of the last replica for which the weight
##          changes. Weights of replicas with indices larger than this index
##          are set to sim.weight_schedule.final

sim.replica.weight_schedule.initial = 1.0
sim.replica.weight_schedule.final   = 0.05
sim.replica.weight_schedule.slope   = 0.001
sim.replica.weight_schedule.first   = 1
sim.replica.weight_schedule.last    = sim.replica.n_replicas

## Generalised temperature for physical prior distribution
##
## initial: (float) Initial prior weight; default: 1.0, value must lie in
##          interval [1., 1.2]
## final:   (float) Final prior weight; default: 1.1
## slope:   (float) Slope of power law used to calculate intermediate weights.
## first:   (integer) Index of the first replica for which the weight begins
##          to change. Weights of replicas with indices smaller than
##          sim.prior_schedule.first are set to sim.prior_schedule.initial
## last:    (integer) Index of the last replica for which the weight changes.
##          Weights of replicas with indices larger than
##          sim.prior_schedule.last are set to sim.prior_schedule.final

sim.replica.prior_schedule.initial = 1.0
sim.replica.prior_schedule.final   = 1.1
sim.replica.prior_schedule.slope   = 0.1
sim.replica.prior_schedule.first   = 1
sim.replica.prior_schedule.last    = sim.replica.n_replicas

###############################################################################
##                                                                           ##
## Settings for analyses / report                                            ##
##                                                                           ##
###############################################################################

## burnin:            (integer) Number of most recent samples to be analysed.
##                    Samples of the initial convergence phase are discarded.
##
## Report:
##
## auto:              (boolean) If set to True, ISD creates a report each time
##                    the results are written to disk (cf. setting
##                    replica.save_interval). Default: False
## keep_sources:      (boolean) Set to True to keep Latex sources (in directory
##                    [WORKING_PATH/analysis/data]. Default: False
##
## pdb_viewer:        (string) Program used to display structures stored in PDB
##                    format. The viewer is called when using the ISD command
##                    "show", with a PDB file as argument
## pdf_latex:         (string) Binary of pdf-latex (default: 'pdflatex').
## eps_to_pdf:        (string) Binary of epstopdf (default: 'epstopdf').
## eps_to_eps:        (string) Binary of eps2eps. Required on some systems,
##                    in order to get the layout of the figures right.
##                    Default: "" (not used).
## max_samples:       (integer) For internal use; should not be modified.
##
## WhatIf:
##
## binary:            (string) Location of the program WhatIf
## use_whatcheck:     (boolean) Set to True if you wish to use WhatCheck
##                    instead of WhatIf. If so, you need to specify the name of
##                    the binary of WhatCheck. Default: False
## show_traces:       (boolean) Set to True in order to generate figures
##                    showing the scores on a per samples basis. Default: False
## enabled:           (boolean) Set to False to disable calculation of WhatIf
##                    quality scores (default: False).
##
## Procheck:
##
## binary:            (string) Location of the program Procheck
## show_traces:       (boolean) Set to True in order to generate figures
##                    showing the scores on a per samples basis. Default: False
## enabled:           (boolean) Set to False to disable calculation of Procheck
##                    quality scores (default: False).
##
## DSSP:
##
## binary:            (string) Location of the program DSSP.
## enabled:           (boolean) Set to False to disable calculation of DSSP
##                    secondary structure assignments (default: False).

sim.analysis.burnin                 = 1000

sim.analysis.report.auto            = False
sim.analysis.report.keep_sources    = False

sim.analysis.pdb_viewer             = 'rasmol'
sim.analysis.gnuplot                = 'gnuplot'
sim.analysis.pdf_latex              = 'pdflatex'
sim.analysis.eps_to_pdf             = 'epstopdf'
sim.analysis.eps_to_eps             = ''
sim.analysis.max_samples            = 1000

sim.analysis.whatif.binary          = 'DO_WHATCHECK.COM'
sim.analysis.whatif.use_whatcheck   = True
sim.analysis.whatif.show_traces     = True
sim.analysis.whatif.enabled         = True

sim.analysis.procheck.binary        = 'procheck'
sim.analysis.procheck.show_traces   = True
sim.analysis.procheck.enabled       = True

sim.analysis.dssp.binary            = 'dssp'
sim.analysis.dssp.enabled           = True

###############################################################################
##                                                                           ##
## Experimental data                                                         ##
##                                                                           ##
###############################################################################
##                                                                           ##
## NOE data                                                                  ##
##                                                                           ##
## Theory: ISPA, i.e. NOE intensity is calculated as inverse six-th power of ##
##         interatomic distances. Ambiguous NOEs are calculated via          ##
##         summation of partial volumes.                                     ##
##                                                                           ##
## Theory parameters: scale of NOEs                                          ##
##                                                                           ##
## Error model: Lognormal distribution with unknown error.                   ##
##                                                                           ##
## In order to maximize the quality of your structure, make sure you add     ##
## your spectra separately via individual data models (see note below).      ##
## That is, avoid merging your peaks into one list. Also, do not remove      ##
## intra-residual NOEs.                                                      ##
##                                                                           ##
## Note: In order to incorporate n NOE data sets, simply duplicate the       ##
##       the block below n times.                                            ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = NOESY
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line.
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## scale.update:  (boolean) Set to True (default) to estimate the scale of the
##                NOEs
## error.initial: Initial value of the error (relative error on distance
##                scale in percent)
## error.update:  (boolean) Set to True (default) to estimate the error of your
##                data

specs.theory.scale.update       = True
specs.error_model.error.initial = 20.
specs.error_model.error.update  = True

##
## Uncomment the following line to add your data to the simulation.
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Distance data                                                             ##
##                                                                           ##
## Theory: None                                                              ##
##                                                                           ##
## Theory parameters: Scale                                                  ##
##                                                                           ##
## Error model: Lognormal distribution with unknown error.                   ##
##                                                                           ##
## Note: In order to incorporate n distance data sets, simply duplicate the  ##
##       the block below n times.                                            ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = DISTANCE
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line.
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## scale.update:  (boolean) Set to True to estimate the scale of distances;
##                default: False
## error.initial: Initial value of the error in percent.
## error.update:  (boolean) Set to True (default) to estimate the error of your
##                data

specs.theory.scale.update       = False
specs.error_model.error.initial = 20.
specs.error_model.error.update  = True

##
## Uncomment the following line to add your data to the simulation
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Scalar couplings                                                          ##
##                                                                           ##
## Theory: Karplus-curve                                                     ##
##                                                                           ##
## Theory parameters: coefficients A, B, C of Karplus curve                  ##
##                                                                           ##
## Error model: Normal distribution with unknown error.                      ##
##                                                                           ##
## Note: In order to incorporate n J coupling data sets, simply duplicate    ##
##       the block below n times.                                            ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = JCOUPLING
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## A, B, C:              (float) Values of the coefficients for the Karplus
##                       relationship
## karplus_curve.update: (boolean) Set to True (default) to estimate
##                       coefficients A, B, C.
## error.initial:        Initial value of the error.
## error.update:         (boolean) Set to True (default) to estimate the error
##                       of your data

specs.theory.karplus_curve.A      = 2.
specs.theory.karplus_curve.B      = -10.
specs.theory.karplus_curve.C      = 1.
specs.theory.karplus_curve.update = True

specs.error_model.error.initial   = 1.
specs.error_model.error.update    = True

##
## Uncomment the following line to add your data to the simulation.
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Dipolar couplings                                                         ##
##                                                                           ##
## Theory: Saupe tensor                                                      ##
##                                                                           ##
## Theory parameters: tensor elements, s1, s2, s3, s4, s5                    ##
##                                                                           ##
## Error model: Normal distribution with unknown error.                      ##
##                                                                           ##
## Note: In order to incorporate n dipolar coupling data sets, simply        ##
## duplicate the the block below n times.                                    ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = RDC
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## s1, s2, s3, s4, s5:   (float) Values of the tensor elements
## saupe_tensor.update:  (boolean) Set to True (default) to estimate the tensor
##                       elements
## error.initial:        Initial value of the error.
## error.update:         (boolean) Set to True (default) to estimate the error
##                       of your data

specs.theory.saupe_tensor.s1     = 1.
specs.theory.saupe_tensor.s2     = 1.
specs.theory.saupe_tensor.s3     = 0.
specs.theory.saupe_tensor.s4     = 0.
specs.theory.saupe_tensor.s5     = 0.
specs.theory.saupe_tensor.update = True

specs.error_model.error.initial  = 1.
specs.error_model.error.update   = True

##
## Uncomment the following line to add your data to the simulation.
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Dihedral angles                                                           ##
##                                                                           ##
## Theory: None                                                              ##
##                                                                           ##
## Theory parameters: None                                                   ##
##                                                                           ##
## Error model: Von-Mises distribution with unknown error.                   ##
##                                                                           ##
## Note: In order to incorporate n dihedral angle data sets, simply          ##
##       duplicate the block below n times.                                  ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, TALOS, or CCPN;
##                default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = DIHEDRAL
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## error.initial: Initial value of the error.
## error.update:  (boolean) Set to True (default) to estimate the error of your
##                data.

specs.error_model.error.update  = True

##
## Uncomment the following line to add your data to the simulation.
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Hydrogen bonds                                                            ##
##                                                                           ##
## Theory: incorporated via distances.                                       ##
##                                                                           ##
## Theory parameters: None                                                   ##
##                                                                           ##
## Error model: Lognormal distribution with unknown error                    ##
##                                                                           ##
## Note: In order to incorporate n hydrogen bond data sets, simply duplicate ##
##       the the block below n times.                                        ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = HBOND
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line.
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## error.initial: Initial value of the error.
## error.update:  (boolean) Set to True (default) to estimate the error of your
##                data

specs.error_model.error.initial = 1.
specs.error_model.error.update  = True

##
## Uncomment the following line to add your data to the simulation
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Disulfide bridges                                                         ##
##                                                                           ##
## Theory: incorporated via distances.                                       ##
##                                                                           ##
## Theory parameters: None                                                   ##
##                                                                           ##
## Error model: Lognormal distribution with unknown error.                   ##
##                                                                           ##
## Note: In order to incorporate n disulfide bridge data sets, simply        ##
##       duplicate the the block below n times.                              ##
##                                                                           ##
###############################################################################

##
## BLOCK START
##
## data_filename: (string) Name of file that contains the data.
## data_format:   (const) Data format; XML, TBL, or CCPN; default: XML
## data_key:      (string) Key to identify data set.
## data_name:     (string or None) Unique name for data set. If set to None,
##                data_key is used as name. Must not contain white space.

data_type     = DISULFIDE
data_filename = ''
data_format   = XML
data_key      = ''
data_name     = None

## Please do not modify this line.
specs = setup_data(data_type, data_filename, data_name, data_key, data_format)

## distance:      (float) SG-SG Distance in Angstrom; default: 2.02
## error.update:  (boolean) Disulfide bridges are enforced. Set to True to
##                estimate the error of your data; default: False

specs.distance                  = 2.02
specs.error_model.error.update  = False

## error.initial: (float) Strength of disulfide bridge; default: 5.0

specs.error_model.error.initial = 5.

##
## Uncomment the following line to add your data to the simulation
##
## sim_data.append(specs)
##
## BLOCK END
##

###############################################################################
##                                                                           ##
## Miscellaneous settings                                                    ##
##                                                                           ##
###############################################################################

## use_xterm: (boolean) For debugging only; if set to True, all replicas are
##            displayed in individual xterm windows. This requires X
##            forwarding to be activated.

sim.use_xterm = False

###############################################################################
##                                                                           ##
## The instance of a simulation can be accessed directly by modifying the    ##
## function below.                                                           ##
##                                                                           ##
## Note: You should know what you are doing here ...                         ##
##                                                                           ##
###############################################################################

def modify_simulation(simulation):
    pass

###############################################################################
##                                                                           ##
## End of setup                                                              ##
##                                                                           ##
## Please do not modify the lines below.                                     ##
##                                                                           ##
###############################################################################

## CCPN: sim.finalise will be called from IsdFrame
## sim.finalize()

if __name__ == '__main__':

    manager = SetupManager(sim, modify_simulation)

    simulation = manager.create_simulation()

    chain = simulation.posterior.get_polymer()

    print 'Starting calculation...'

    sampler = manager.create_sampler()

    initial_states = manager.prepare_initial_states(sampler)

    print 'Calculation started.'

    sampler.generate_sequence(sim.replica.n_samples, initial_states)


