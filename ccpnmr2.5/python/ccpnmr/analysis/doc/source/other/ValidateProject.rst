==================
Project : Validate
==================

**Check the CCPN project for data consistency**

This option performs an analysis of the data stored in the current CCPN project and
checks its consistency with respect to the description of the CCPN data model.

Although the API to the CCPN data model continuously performs checks whenever data is
manipulated in Analysis, it is possible for the project information to become
inconsistent within itself if the program or computer suffered an error at a critical
time. This may occur for example if the program crashes in the middle of a save or if an
error ocurs in Python code that bypasses the main API routines.

The user may activate the project validation function to give peace of mind that a
project will be saved properly, or has been loaded properly, after an error interrupted
the program.

It should be noted that this check only tests the sanity of the project with regards to
the definition of the objects, links and attributes within the CCPN data model. It does
not check the scientific content of the project for errors. For identifying some
scientific errors in NMR data analysis the `Quality Reoprts`_ option may be used. For
identifying errors that relate to the generation of macromolcular structures the CING_
option can be used to validate via the iCing server.

.. _`Quality Reoprts`: ../popups/ViewQualityReportsPopup.html
.. _CING: ../popups/CingPopup.html
