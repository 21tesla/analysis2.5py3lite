"""Shared fixtures for CCPNMR core test suite."""
import os
import sys

# Ensure the ccpnmr2.5/python path is on sys.path (belt-and-braces for pytest)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
