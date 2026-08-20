CCPN Data Model API 1.0 release notes

A short tour of the important parts of the release, by directory:

ccpnmr/                            top level directory

  INSTALL                          installation notes
  README                           release notes (this file)
  installCode.py                   Python script to install code

  bin/                             executables

  doc/                             top level documentation

  ccpnmr1.0/                       CcpNmr code

    python/                        Python code

      doc/                         data model API documentation

      memops/                      non-NMR code

        universal/                 utility code

        general/                   data model utility code

        api/                       basic data model code
        format/
        metamodel/
        model/
        xml/

      ccp/                         general NMR code

        general/                   data model utility code
        util/

        api/                       basic data model code
        model/
        xml/


    data/                          some data (e.g. molecular definitions)

  python2.2/                       third-party Python code

  tcl8.3/                          third-party Tcl code

  tk8.3/                           third-party Tk code

  mesa4.0/                         third-party OpenGL (Mesa) code
