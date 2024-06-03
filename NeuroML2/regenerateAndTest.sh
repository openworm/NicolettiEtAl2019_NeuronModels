#!/bin/bash
set -ex

black *.py

python GenerateNeuroML.py -jnml

python TestXPP.py

python GenerateExamples.py -jnml

omv all -V 

