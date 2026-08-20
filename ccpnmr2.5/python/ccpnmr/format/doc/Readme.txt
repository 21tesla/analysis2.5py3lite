CcpNmr FormatConverter 2.1 release notes

Assuming successful installation, to start FormatConverter type:

  formatConverter [ project_file ]

where the project_file is an optional argument specifying a
previously saved project.

A short tour of the important parts of the release, by directory:

ccpnmr/                            top level directory

  INSTALL                          installation notes
  README                           release notes (this file)
  installCode.py                   Python script to install code

  bin/                             executables

  doc/                             top level documentation

  ccpnmr2.1/                       CcpNmr code

    python/                        Python code

      doc/                         data model API documentation

      memops/                      non-NMR code

        universal/                 utility code

        general/                   data model utility code

        gui/                       Tkinter widget code

        editor/                    data model Tkinter widget code

        math/                      math-related code

        api/                       basic data model code
        format/
        metamodel/
        xml/

      ccp/                         general NMR code

        general/                   data model utility code
        util/

        math/                      math-related code

        format/                    NMR format conversion code

        api/                       basic data model code
        xml/

      ccpnmr/                      CcpNmr-specific code

        format/                    CcpNmr format conversion code

        api/                       basic data model code
        xml/

      pdbe/                        PDBE (deposition) code

      molsim/                      molecular simulation code

        api/                       basic data model code
        xml/

    data/                          some data (e.g. molecular definitions)

  python2.6/                       third-party Python code (optional)

  tcl8.5/                          third-party Tcl code (optional)

  tk8.5/                           third-party Tk code (optional)
