#!/bin/bash
set -ex

black *.py

python GenerateNeuroML.py -jnml

omv all -V 

