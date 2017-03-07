#!/bin/bash
# too lazy to make a Makefile

nvcc --ptxas-options=-v -O3 -o fault_sim fault_sim.cu

nvcc --ptxas-options=-v -O3 --ptx --source-in-ptx -o fault_sim.ptx fault_sim.cu

