#!/bin/bash
set -ex

# Format the code
black *.py

# This script loads in the XPP using the parser in pyNeuroML, and converts the channels to NeuroML
# It will also run a test simulation of the cells and save the results
# NOTE: make sure to update to a recent version of pyNeuroML!
python GenerateNeuroML.py -jnml -short


# This will load the original XPP using pyNeuroML, run it, save the results and plot them against the 
# results generated by the NeuroML in GenerateNeuroML.py
python TestXPP.py -short


