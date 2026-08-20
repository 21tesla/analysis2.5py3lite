""" adapter to import the best ElementTree implementation
dependent on Python version
"""

try:
  from lxml import etree as ElementTree
except ImportError:
  try:
    import xml.etree.cElementTree as ElementTree # in python >=2.5
  except ImportError:
    try:
      from  xml.etree import ElementTree # in python >=2.5
    except ImportError:
      from elementtree import ElementTree # effbot's pure Python module

try:
  from lxml import ElementInclude
except ImportError:
  try:
    from  xml.etree import ElementInclude # in python >=2.5
  except ImportError:
    from elementtree import ElementInclude # effbot's pure Python module
