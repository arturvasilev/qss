#!/bin/bash
cd ~/qss
mpic++ -std=c++11 main.cpp graph.cpp -lnlopt -lm -o quantum_scheme_solver
#./quantum_scheme_solver