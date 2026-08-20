Making Bug Reports
==================

Only the latest software version is actively supported, but old CCPN data files
are supported.

Consider updating Analysis to the latest version before reporting an issue that
may have been fixed. Remember to include any recent updates or patches: Analysis
is sometimes updated on a daily basis.

Bugs will often be fixed promptly, but if the resolution is more involved they
may be put on the to-do list.

If a bug seems important, but seems to be neglected at the expense of other
issues then it has probably slipped though the cracks and the developers should
asked what the status is.

In essence, state what you did, with enough information for the developers to
understand and reproduce the problem. The better the information the more
likely the problem is quick to fix.

Information to Include
----------------------

Give the name of the popup window that shows any inconsistencies or was being
used when the error occurred.

Report any errors that appear on the command line, even if they occurred some
time before you noticed anything was wrong. If you are not sure whether
something is a bug, you can always ask.

Be precise, don't say:

  *"The shift list doesn't show my molecule assignments properly"*

when you really mean:

  *"The Chart:Chemical Shifts Table, alone, doesn't show backbone 15N assignments
  for some residues in chain B."*

For NMR concepts, try to report using the same sort of language that appears
in the program, e.g. distinguish between "experiment" and "spectrum" or
"molecule" and "chain" when Analysis does.

If it helps, illustrate the problem by including screen shots as attachments,
and consider selecting the relevant rows in any tables.

Broken or illustrative CCPN projects should be sent as a tar-gzipped archive,
without spectrum data files.

Reporting Style
---------------

Try not to make assumptions about whether something is a bug or design decision.
Sometimes a design decision can appear like a bug from some points of view, and
some bugs look like features. Often it is best to simply ask what was intended.

Try not to prescribe *how* something should be fixed unless you are proposing to
fix something yourself. Most problems are best described in terms of the
functional requirements.

Suggestions for functionality and changes may be made at any time, but repeating
yourself should the developers or the user community disagree with the
suggestions is mostly futile.

Dwelling on the hardship a problem has caused may delay the response rather than
hasten it.
