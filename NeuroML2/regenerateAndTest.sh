#!/bin/bash
set -ex

black *.py

python GenerateNeuroML.py -jnml

python TestXPP.py

omv all -V 

