#!/bin/bash

INPUT_GEO_FILE="domain.geo"
INITIAL_MSH_FILE="domain.msh"

gmsh $INPUT_GEO_FILE -2 -o $INITIAL_MSH_FILE
