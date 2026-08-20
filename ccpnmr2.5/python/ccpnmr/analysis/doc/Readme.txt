CcpNmr Analysis 2.4 release notes

Assuming successful installation, if ccpnmr/bin is on your PATH, then
to start Analysis type:

  analysis [ project_file ]

where the project_file is an optional argument specifying a
previously saved project.

A short tour of the important parts of the release, by directory:

ccpnmr/                            top level directory

  INSTALL                          installation notes
  README                           release notes (this file)
  installCode.py                   Python script to install code

  bin/                             scripts ("executables")

  doc/                             top level documentation

  ccpnmr2.5/                       CcpNmr code (this folder name is subject to change)

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

        gui/                       GUI code that also involves data model

        math/                      math-related code

        format/                    NMR format conversion code

        api/                       basic data model code
        xml/

      ccpnmr/                      CcpNmr-specific code

        analysis/                  CcpNmr Analysis code

        clouds/                    CcpNmr Clouds code

        nexus/                     CcpNmr Nexus code

        eci/                       CcpNmr ECI code

        format/                    CcpNmr format conversion code

        workflow/                  CcpNmr workflow related code

        api/                       basic data model code
        xml/

      pdbe/                        PDBE (deposition) code

      molsim/                      molecular simulation code

        api/                       basic data model code
        xml/

      model/                       data model

      cambridge/                   third-party NMR-related code
      extendNmr/
      gothenburg/
      gottingen/
      nijmegen/
      paris/
      regensburg/
      utrecht/

    c/                             C code

    data/                          some data (e.g. chemical component definitions)

The remaining directories are optional (so come from separate installs):

  python2.6/                       third-party Python code

  numpy1.4/                        third-party NumPy code

  tcl8.5/                          third-party Tcl code

  tk8.5/                           third-party Tk code

  tix8.4/                          third-party Tix code

  mesa6.0/                         third-party OpenGL (Mesa) code

