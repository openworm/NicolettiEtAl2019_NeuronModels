#!/bin/bash
set -ex

python GenerateNeuroML.py -jnml

omv all -V 

