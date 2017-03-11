#!/bin/bash

OPTS="-O3 --gpu-architecture=compute_50 --gpu-code=compute_50,sm_50,sm_52"
nvcc $OPTS -DUSE_TRAX=1 -DUSE_HAZARDS=1 -o fault_sim_trax fault_sim.cu
nvcc $OPTS -DUSE_TRAX=1 -o fault_sim_traxnh fault_sim.cu
nvcc $OPTS -DUSE_TRAX=0 -o fault_sim_tf fault_sim.cu

# This is for debugging?
#nvcc --ptxas-options=-v -O3 --ptx --source-in-ptx -o fault_sim.ptx fault_sim.cu

