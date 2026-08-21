'''
Created on Jun 24, 2011

@author: jd
'''

try:
    import nose
except ModuleNotFoundError:
    nose = None

if nose.run():
    print("Nose ran fine")
else:
    print("ERROR: Nose failed")
