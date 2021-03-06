/*

   GPU fault simulator for TRAX fault model
   Matthew Beckler
   18-645 How to Write Fast Code
   Carnegie Mellon University
   Fall 2012

   Last Updated: March 22, 2017
   Copyright (c) 2013-2017, Matthew Beckler

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


   Usage: ./fault_sim basename
   Looks for the files basename/basename.easy and basename/basename.tests and TODO
   Stores dictionary output in dictionary.cuda

   For more information about this GPU fault simulator:

   M. Beckler and R. D. S. Blanton, "GPU-accelerated fault dictionary generation for the TRAX fault model," 2017 International Test Conference in Asia (ITC-Asia), Taipei, 2017, pp. 34-39. - doi: 10.1109/ITC-ASIA.2017.8097107 - https://doi.org/10.1109/ITC-ASIA.2017.8097107

   M. Beckler and R. D. Blanton, "Fault simulation acceleration for TRAX dictionary construction using GPUs," 2017 IEEE International Test Conference (ITC), Fort Worth, TX, 2017, pp. 1-9. - doi: 10.1109/TEST.2017.8242078 - https://doi.org/10.1109/TEST.2017.8242078

   For more information about the TRAX fault model:

   M. Beckler and R. D. Blanton, "On-Chip Diagnosis of Generalized Delay Failures using Compact Fault Dictionaries," in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. PP, no. 99, pp. 1-1. - doi: 10.1109/TCAD.2018.2803621 - https://doi.org/10.1109/TCAD.2018.2803621

   Beckler, M. and Blanton, R.D., "On-Chip Diagnosis for Early-Life and Wear-Out Failures," IEEE International Test Conference, Nov. 2012. - https://doi.org/10.1109/TEST.2012.6401580


   TODO list of things to improve:

   * Figure out what's going on with dev_gates, didn't we put this in constants memory?
        Nope, only the gate eval LuT is in constants memory. There's only 64k of constants mem, and we need 128 bytes for the gate eval LuT. Each Gate takes up 12 bytes (and we could only reduce that by limiting us to 2^16 = 65536 nets which might not be ok) so we can have a maximum of 5450 gates, which is not enough. Unless we could use multiple kernels for each chunk of 5450 gates? Probably not worth it.

   * Maybe use shared memory for the netlist information?
        Probably not enough room for the larger circuits we want to handle

   * Store dictionary in memory/disk as binary format?

   * Some way to do multiple gate evals per thread? Maybe a packed data structure? Maybe machine word-width algebraic operations like in that one DAC submission?

   * As mentioned above, some data storage techniques had to be migrated to a compacted storage technique, where multiple 1 or 2 bit values are stored in a single byte. This reduces the total storage required, but makes each access more complex, usually involving one or more bit shift and bitwise masking operations. It would be good to analyze the actual access patterns for these memories, and determine if functions such as "extract four consecutive values" or "extract both v1 and v2 for a given net" would be useful and faster than individual accesses. Additionally, for smaller circuits or on systems with more available memory, it would be useful to be able to automatically detect the memory size and switch between data storage formats based on the system details.

   * There are a number of small optimizations we would like to investigate with regards to array-of-struct vs struct-of-arrays, especially in our gate and test memory storage.

   Update log:
   * Ignore the skip-ahead optimization for Kernel 3, and the N invocations of Kernel 3, by doing one giant Kernel 3 invocation.
   * Support tracking activity of each net (and therefore PMOS input) over the entire test set - Mar 2017
   * Added support for disabling TRAX hazard activation - Mar 2017
   * Added support for TF faults in addition to just TRAX faults - Aug 2013

*/

#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>

/* Here is the "easy" to parse file format for the circuit netlist:
    NUMINPUTS 5
    INPUT 8
    INPUT 9
    INPUT 10
    INPUT 11
    INPUT 12
    NUMOUTPUTS 4
    OUTPUT 6
    OUTPUT 7
    OUTPUT 1
    OUTPUT 0
    NUMGATES 8
    1 0 9 9
    0 1 8 8
    2 2 8 10
    3 3 10 11
    4 4 9 3
    5 5 3 12
    8 6 2 4
    9 7 4 5

    Nets are integers.
    Gates can only have two inputs (buf/inv should have the same net listed as both inputs, the second will be ignored, but won't break the parser).
    GATES ABSOLUTELY MUST BE TOPOLOGICALLY SORTED based on their level within the circuit.
    Gates' output net ids must start at 0.
    There are probably other implicit assumptions that I have forgotten about. Caveat emptor
*/

#define uchar unsigned char
#define uint unsigned int

// This is the definition of all eight supported gate types. Ordering is to make the GPU code more efficient (and/nand/or/nor are commonly grouped together, same with xor/xnor)
#define TYPE_AND    (0)
#define TYPE_NAND   (1)
#define TYPE_OR     (2)
#define TYPE_NOR    (3)
#define TYPE_BUF    (4)
#define TYPE_INV    (5)
#define TYPE_XOR    (6)
#define TYPE_XNOR   (7)

// Since we have a four-valued logic, we need two bits to represent each one.
// I was originally trying to be creative about these assignments, to make it so that we could do gate evaluations in parallel using bit-operation instructions.
// Which is how I decided to make 0 = 00, 1 = 11, and X and H be the other two.
// However, since we haven't implemented parallel gate evaluations, this is kind of a moot point. However these values are indeed being used.
#define LOGIC_0 (0)
#define LOGIC_X (1)
#define LOGIC_H (2)
#define LOGIC_1 (3)

// How many threads per block? Maximum is 512, make sure each is a multiple of 16, multiple of 32 is probably better.
// Doesn't really seem to affect performance much, so something else is probably limiting the speed.
#define FAULTS_PER_BLOCK_KERNEL_1 (512)
#define FAULTS_PER_BLOCK_KERNEL_2 (512)
#define FAULTS_PER_BLOCK_KERNEL_3 (512)

// FYI, pointers cost 8 bytes!

// Structure to represent a gate
// We should pack in1 and in2 together in the same int if we need to save memory
typedef struct
{
    uchar type;
    uchar is_output;
    uint in1;
    uint in2;
} Gate; // 12 bytes due to nice packing

// Structure to represent a fault
typedef struct
{
    uint net;
    uint polarity; // 0 = falling, 1 = rising
} Fault;


// Used to evaluate how long (wall-clock) the three kernels take to run. Taken from:
// http://stackoverflow.com/questions/1468596/c-programming-calculate-elapsed-time-in-milliseconds-unix
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}
// Accumulates (end - start) into result
void timeval_accumulate_diff(struct timeval *result, struct timeval *start, struct timeval *end)
{
    long int diff_usec = (end->tv_usec + 1000000 * end->tv_sec) - (start->tv_usec + 1000000 * start->tv_sec);
    diff_usec += result->tv_usec;
    result->tv_sec += diff_usec / 1000000;
    result->tv_usec = diff_usec % 1000000;
}

// Prints a nice number of bytes (KB, MB, GB, etc) to the provided buffer.
void pretty_bytes(char* buf, uint bytes)
{
    const char* suffixes[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};
    uint order = 0;
    double count = bytes;
    while ((count >= 1024) && (order < sizeof(suffixes))) {
        order++;
        count /= 1024;
    }
    if ((count - floor(count)) == 0.0) {
        sprintf(buf, "%d %s", (int)count, suffixes[order]);
    } else {
        sprintf(buf, "%.1f %s", count, suffixes[order]);
    }
}


char printable_logic_value(char value)
{
    switch (value)
    {
        case LOGIC_0: return '0';
        case LOGIC_X: return 'X';
        case LOGIC_H: return 'H';
        case LOGIC_1: return '1';
    }
    return '?';
}


// Circuit states are stored as pairs of values (v1, v2) for each net in the circuit.
void print_state(uchar* state, uint num_state_values)
{
    for (uint i = 0; i < num_state_values; i++)
    {
        uchar s = state[i];
        printf("%c", printable_logic_value(s));
        if (i % 2 == 1)
            printf(" ");
    }
    printf("\n");
}
void print_state_raw(uchar* state, uint num_state_values)
{
    for (uint i = 0; i < num_state_values; i++)
    {
        printf("%X", state[i]);
    }
    printf("\n");
}

// For some data like the tests values we pack 8 1-bit values into each byte.
// These are the accessor functions to make it easy to access that data.
__host__ __device__ void BIT_SET_UCHAR(uchar* array, uint which, uchar value)
{
    uchar* b = array + (which / 8);
    uchar shift = which % 8;
    if (value)
        *b |= (1 << shift);
    else
        *b &= ~(1 << shift);
}
__host__ __device__ uchar BIT_GET_UCHAR(uchar* array, uint which)
{
    uchar* b = array + (which / 8);
    uchar shift = which % 8;
    if (*b & (1 << shift))
        return 1;
    else
        return 0;
}

// This is the new gate evaluation lookup table.
// Simply concatenate the gate type with the four bits of input values.
__constant__ uchar gate_eval_lut[128] = {
    // gate 000 (and)
    LOGIC_0, // 0 0 -> 0
    LOGIC_0, // 0 X -> 0
    LOGIC_0, // 0 H -> 0
    LOGIC_0, // 0 1 -> 0
    LOGIC_0, // X 0 -> 0
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_X, // X 1 -> X
    LOGIC_0, // H 0 -> 0
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H
    LOGIC_0, // 1 0 -> 0
    LOGIC_X, // 1 X -> X
    LOGIC_H, // 1 H -> H
    LOGIC_1, // 1 1 -> 1

    // gate 001 (nand)
    LOGIC_1, // 0 0 -> 1
    LOGIC_1, // 0 X -> 1
    LOGIC_1, // 0 H -> 1
    LOGIC_1, // 0 1 -> 1
    LOGIC_1, // X 0 -> 1
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_X, // X 1 -> X
    LOGIC_1, // H 0 -> 1
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H
    LOGIC_1, // 1 0 -> 1
    LOGIC_X, // 1 X -> X
    LOGIC_H, // 1 H -> H
    LOGIC_0, // 1 1 -> 0

    // gate 010 (or)
    LOGIC_0, // 0 0 -> 0
    LOGIC_X, // 0 X -> X
    LOGIC_H, // 0 H -> H
    LOGIC_1, // 0 1 -> 1
    LOGIC_X, // X 0 -> X
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_1, // X 1 -> 1
    LOGIC_H, // H 0 -> H
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_1, // H 1 -> 1
    LOGIC_1, // 1 0 -> 1
    LOGIC_1, // 1 X -> 1
    LOGIC_1, // 1 H -> 1
    LOGIC_1, // 1 1 -> 1

    // gate 011 (nor)
    LOGIC_1, // 0 0 -> 1
    LOGIC_X, // 0 X -> X
    LOGIC_H, // 0 H -> H
    LOGIC_0, // 0 1 -> 0
    LOGIC_X, // X 0 -> X
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_0, // X 1 -> 0
    LOGIC_H, // H 0 -> H
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_0, // H 1 -> 0
    LOGIC_0, // 1 0 -> 0
    LOGIC_0, // 1 X -> 0
    LOGIC_0, // 1 H -> 0
    LOGIC_0, // 1 1 -> 0

    // gate 100 (buffer)
    LOGIC_0, // 0 0 -> 0
    LOGIC_0, // 0 X -> 0 - not possible, for BUF/INV must have same values
    LOGIC_0, // 0 H -> 0 - not possible, for BUF/INV must have same values
    LOGIC_0, // 0 1 -> 0 - not possible, for BUF/INV must have same values
    LOGIC_X, // X 0 -> X - not possible, for BUF/INV must have same values
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X - not possible, for BUF/INV must have same values
    LOGIC_X, // X 1 -> X - not possible, for BUF/INV must have same values
    LOGIC_H, // H 0 -> H - not possible, for BUF/INV must have same values
    LOGIC_H, // H X -> H - not possible, for BUF/INV must have same values
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H - not possible, for BUF/INV must have same values
    LOGIC_1, // 1 0 -> H - not possible, for BUF/INV must have same values
    LOGIC_1, // 1 X -> H - not possible, for BUF/INV must have same values
    LOGIC_1, // 1 H -> H - not possible, for BUF/INV must have same values
    LOGIC_1, // 1 1 -> H

    // gate 001 (inverter)
    LOGIC_1, // 0 0 -> 1
    LOGIC_1, // 0 X -> 1 - not possible, for BUF/INV must have same values
    LOGIC_1, // 0 H -> 1 - not possible, for BUF/INV must have same values
    LOGIC_1, // 0 1 -> 1 - not possible, for BUF/INV must have same values
    LOGIC_X, // X 0 -> X - not possible, for BUF/INV must have same values
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X - not possible, for BUF/INV must have same values
    LOGIC_X, // X 1 -> X - not possible, for BUF/INV must have same values
    LOGIC_H, // H 0 -> H - not possible, for BUF/INV must have same values
    LOGIC_H, // H X -> H - not possible, for BUF/INV must have same values
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H - not possible, for BUF/INV must have same values
    LOGIC_0, // 1 0 -> 0 - not possible, for BUF/INV must have same values
    LOGIC_0, // 1 X -> 0 - not possible, for BUF/INV must have same values
    LOGIC_0, // 1 H -> 0 - not possible, for BUF/INV must have same values
    LOGIC_0, // 1 1 -> 0

    // gate 110 (xor)
    LOGIC_0, // 0 0 -> 0
    LOGIC_X, // 0 X -> X
    LOGIC_H, // 0 H -> H
    LOGIC_1, // 0 1 -> 1
    LOGIC_X, // X 0 -> X
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_X, // X 1 -> X
    LOGIC_H, // H 0 -> H
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H
    LOGIC_1, // 1 0 -> 1
    LOGIC_X, // 1 X -> X
    LOGIC_H, // 1 H -> H
    LOGIC_0, // 1 1 -> 0

    // gate 111 (xnor)
    LOGIC_1, // 0 0 -> 1
    LOGIC_X, // 0 X -> X
    LOGIC_H, // 0 H -> H
    LOGIC_0, // 0 1 -> 0
    LOGIC_X, // X 0 -> X
    LOGIC_X, // X X -> X
    LOGIC_X, // X H -> X
    LOGIC_X, // X 1 -> X
    LOGIC_H, // H 0 -> H
    LOGIC_X, // H X -> X
    LOGIC_H, // H H -> H
    LOGIC_H, // H 1 -> H
    LOGIC_0, // 1 0 -> 0
    LOGIC_X, // 1 X -> X
    LOGIC_H, // 1 H -> H
    LOGIC_1, // 1 1 -> 1
};


// This is the core cuda fault simulation function, only called from the two fault sim wrapper functions (kernels 1 and 3).
// It will probably be inlined into each kernel by the compiler, but this is better for future maintenance
// Inputs:
//      gates           The global, unchanging gate structure
//      my_state        Pointer to the state for this simulation (explained above, basically it's just net0_v1, net0_v2, net1_v1, net1_v2, etc)
// Outputs:
//      Updates my_state in place with the correct values after fault simulation
// Returns:
//      The v2 value
#define TRAN0110 ( (LOGIC_0 << 6) | (LOGIC_1 << 4) | (LOGIC_1 << 2) | LOGIC_0 )
#define TRAN1001 ( (LOGIC_1 << 6) | (LOGIC_0 << 4) | (LOGIC_0 << 2) | LOGIC_1 )
#define TRAN01H0 ( (LOGIC_0 << 6) | (LOGIC_1 << 4) | (LOGIC_H << 2) | LOGIC_0 )
#define TRAN100H ( (LOGIC_1 << 6) | (LOGIC_0 << 4) | (LOGIC_0 << 2) | LOGIC_H )
#define TRAN011H ( (LOGIC_0 << 6) | (LOGIC_1 << 4) | (LOGIC_1 << 2) | LOGIC_H )
#define TRAN10H1 ( (LOGIC_1 << 6) | (LOGIC_0 << 4) | (LOGIC_H << 2) | LOGIC_1 )
__device__ uchar cuda_fault_sim_core(const Gate* g, uchar* my_state, uint gate_id)
{
    // could we transpose our state matrix to make the accesses align better? TODO
    // Each thread is accessing the same net at the same time, maybe we could make those accesses be coallesced?
    uchar in1_v1 = 0x7F & my_state[g->in1 * 2];
    uchar in2_v1 = 0x7F & my_state[g->in2 * 2];
    uchar in1_v2 = 0x7F & my_state[g->in1 * 2 + 1];
    uchar in2_v2 = 0x7F & my_state[g->in2 * 2 + 1];

    // Using the new gate evaluation lookup table, we use the gate type and two input values to craft a 7-bit index
    uchar v1 = gate_eval_lut[(g->type << 4) | (in1_v1 << 2) | (in2_v1)];
    uchar v2 = gate_eval_lut[(g->type << 4) | (in1_v2 << 2) | (in2_v2)];

#if (USE_TRAX && USE_HAZARDS)
    // Now we need to detect if v2 should result in a hazard due to this gate's inputs
    uchar merged_values = (in1_v1 << 6) | (in2_v1 << 4) | (in1_v2 << 2) | in2_v2;
    // and/nand = 0/1, or/nor = 2/3
    uchar hazard_and_nand = (g->type <= 1) &&                 (merged_values == TRAN0110 || merged_values == TRAN1001 || merged_values == TRAN01H0 || merged_values == TRAN100H);
    uchar hazard_or_nor   = (g->type <= 3 && g->type >= 2) && (merged_values == TRAN0110 || merged_values == TRAN1001 || merged_values == TRAN011H || merged_values == TRAN10H1);
    uchar hazard_xor_xnor = (g->type >  5) && (v1 == v2 && in1_v1 != in1_v2 && in2_v1 != in2_v2);
    uchar hazard = hazard_and_nand || hazard_or_nor || hazard_xor_xnor;
    v2 = hazard ? LOGIC_H : v2;
#endif

    // Update my_state in place with the (potentially new) values of v1 and v2:
    // Note that unlike the reference implementation, we don't care if the values were updated, since we check all downstream gates regardless.
    my_state[gate_id * 2] = v1;
    my_state[gate_id * 2 + 1] = v2;

    return v2;
}


// This is Kernel 1, the cuda fault-free fault simulation function. Pass in the gates, all the states (one per test), the state size, and the number of gates.
// Updates all the states in place. One thread per test pair, FAULTS_PER_BLOCK_KERNEL_1 threads per block.
__global__ void cuda_fault_free_fault_sim(Gate* gates, uint num_gates, uchar* all_states, uint state_bytes, uint num_tests)
{
    // each thread will go through a complete fault simulation
    // since we are not skipping around, all threads stay in lock-step
    // Key idea: branching is ok as long as all the threads do the same thing
    // This means we can't have any data-dependent branching
    // Since each thread evaluates the same gate at the same time, that kind of branching is ok
    // However we can't have branching based on the actual logic values going through the gate

    uint test_id = (blockIdx.x * FAULTS_PER_BLOCK_KERNEL_1) + threadIdx.x;
    if (test_id < num_tests)
    {
        uchar* my_state = all_states + (state_bytes * test_id);
        // We iterate over all the relevant gates, in topological order (this sorting was already handled, our netlist comes pre-sorted)
        for (uint gate_id = 0; gate_id < num_gates; gate_id++)
        {
            cuda_fault_sim_core(&gates[gate_id], my_state, gate_id);
        }
    }
}

// This is Kernel 2, which determines which tests activate each fault. One thread per fault.
// It stores a 0 or 1 into the matrix fault_activation, which now uses a packed format (one bit for every fault-test pair), each row is a fault, each col is a test. (initialized to all zeros).
__global__ void cuda_check_fault_activations(Gate* gates, uchar* all_states, uint state_bytes, uint num_tests, uchar* fault_activations, uint num_faults)
{
    uint fault_id = blockIdx.x * FAULTS_PER_BLOCK_KERNEL_2 + threadIdx.x;
    if (fault_id < num_faults)
    {
        uint gate_id = fault_id / 2;
        uint rising = (fault_id % 2); // rising when fault_id is odd
        // this thread has to check for fault activation for fault_id for all tests
        for (uint test_id = 0; test_id < num_tests; test_id++)
        {
            uchar* my_state = all_states + test_id * state_bytes;
            uchar v1 = my_state[gate_id * 2];
            uchar v2 = my_state[gate_id * 2 + 1];
            uchar activated = ( ( rising && (v1 == LOGIC_0) && (v2 == LOGIC_1)) || // output transition activation
                                (!rising && (v1 == LOGIC_1) && (v2 == LOGIC_0)) || // ditto
                                (v2 == LOGIC_H) );                                 // hazard-based activation - And actually, if we are not doing TRAX fault sim, then there can be no LOGIC_H values, so this activation condition is harmless, neat!
            // argh, this is also susceptible to the problem of multiple threads writing to the same byte concurrently
            //BIT_SET_UCHAR(fault_activations, fault_id * num_tests + test_id, (activated ? 1 : 0));
            // (TODO - This comment is talking about BIT_SET_UCHAR, right? Not the atomicOr, right?)

            uint index = fault_id * num_tests + test_id;
            atomicOr( ((uint*)fault_activations) + (index / 32), activated << (index % 32));
        }
    }
}


// This is Kernel 3, which does the faulty fault simulation. One grid of blocks per fault, FAULTS_PER_BLOCK_KERNEL_3 threads per block, one thread per activating test.
// NEW PLAN: Update all the states in place, but never copy the faulty states back to the CPU.
// Instead, we now have the dictionary row in memory and we directly write the pass/fail bit into that memory directly, and then copy it back to the CPU and write it to disk.
__global__ void cuda_faulty_fault_sim(Gate* dev_gates, uint num_gates, uchar* dev_faulty_states, uint state_bytes, uint total_activations)
{
    uint offset = blockIdx.x * FAULTS_PER_BLOCK_KERNEL_3 + threadIdx.x;
    if (offset < total_activations)
    {
        uchar* my_state = dev_faulty_states + state_bytes * offset;
        // Removed fault activation since it now happens on the CPU
        uchar test_failed = 0; // Local copy since we'll be writing it many times, copy it to final byte of my_state at the end

        // We iterate over all the relevant gates, in topological order (this sorting was already handled, our netlist comes pre-sorted)
        for (uint gate_id = 0; gate_id < num_gates; gate_id++)
        {
            Gate *g = &dev_gates[gate_id]; // should be the same for all threads in the kernel
            uchar v2 = my_state[gate_id * 2 + 1]; // Existing V2
            if (!(v2 & 0x80)) {
                // 0x80 set = "activated fault site" so we don't re-simulate it and un-activate the fault!
                v2 = cuda_fault_sim_core(g, my_state, gate_id);
            }

        #if USE_TRAX
            test_failed = (g->is_output && ((0x7F & v2) == LOGIC_X)) ? 1 : test_failed; // If an X reaches an output, the test fails
        #else
            test_failed = (g->is_output && (v2 != my_state[gate_id * 2 + 1])) ? 1 : test_failed; // If an output does not match the expected (fault-free) value, the test fails
        #endif
        }
        my_state[state_bytes - 1] = test_failed;
    }
}

void check_cuda_errors(const char *kernel_name)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error in kernel %s: %s\n", kernel_name, cudaGetErrorString(err));
        exit(1);
    }
}

// does an integer division of n/d, rounding up
uint divide_round_up(uint n, uint d)
{
    return (n + (d - 1)) / d;
}

// TODO use compile-time directives for the TRAX vs TF code to improve performance?

// This is our main function.
// It would be stellar to split the parsing code into a separate function in a separate file.
int main(int argc, char* argv[])
{
    char buffer[100];  // for using pretty_bytes(buffer, numbytes);

    // Track how many bytes we malloc on the CPU and the CUDA device
    unsigned long malloc_bytes = 0;
    unsigned long cuda_malloc_bytes = 0;

    // I added some timing instrumentation in the code for different sections
    struct timeval tvStart, tvDoneParsing, tvPreK1, tvPostK1, tvPreK2, tvPostK2, tvPreK3, tvPostK3, tvEnd, tvDiff;
    gettimeofday(&tvStart, NULL);

    printf("This code built with these recognized options:\n"
#if USE_TRAX
           "USE_TRAX "
#   if USE_HAZARDS
           "USE_HAZARDS "
#   else
           "NO_HAZARDS "
#   endif
#else
           "USE_TF "
#endif
           "\n");

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s basename}\n", argv[0]);
        exit(1);
    }

    char* basename = argv[1]; // something like "c432"

    malloc_bytes += (2 * strlen(basename) + 6 + 1);
    char* filename_input = (char*) malloc(2 * strlen(basename) + 6 + 1); // "c432/c432.v", so 2*N + 6 for "/.easy" + 1 for the '\0'
    assert(filename_input != NULL);
    sprintf(filename_input, "%s/%s.easy", basename, basename);
    printf("Netlist filename: '%s'\n", filename_input);

    malloc_bytes += (2 * strlen(basename) + 12 + 1);
    char* filename_tests = (char*) malloc(2 * strlen(basename) + 12 + 1); // "c432/c432.tests.easy", so 2*N + 12 for "/.tests.easy" + 1 for the '\0'
    assert(filename_tests != NULL);
    sprintf(filename_tests, "%s/%s.tests.easy", basename, basename);
    printf("Tests filename: '%s'\n", filename_tests);

    malloc_bytes += (2 * strlen(basename) + 27 + 1);
    char* filename_dictionary = (char*) malloc(2 * strlen(basename) + 27 + 1); // "c432/c432.dictionary.traxnh.pf.cuda", so 2*N + 27 for "/.dictionary.traxnh.pf.cuda" + 1 for the '\0'
    assert(filename_dictionary != NULL);
    sprintf(filename_dictionary,
#if USE_TRAX
#   if USE_HAZARDS
    "%s/%s.dictionary.trax.pf.cuda",
#   else
    "%s/%s.dictionary.traxnh.pf.cuda",
#   endif
#else
    "%s/%s.dictionary.tf.pf.cuda",
#endif
    basename, basename);
    printf("Dictionary filename: '%s'\n", filename_dictionary);

    malloc_bytes += (2 * strlen(basename) + 13 + 1);
    char* filename_faults = (char*) malloc(2 * strlen(basename) + 13 + 1); // "c432/c432.faults.gpu", so 2*N + 13 for "/.faults.gpu" + 1 for the '\0'
    assert(filename_faults != NULL);
    sprintf(filename_faults, "%s/%s.faults.gpu", basename, basename);
    printf("Faults filename: '%s'\n", filename_faults);

    malloc_bytes += (2 * strlen(basename) + 7 + 1);
    char* filename_usage = (char*) malloc(2 * strlen(basename) + 7 + 1); // "c432/c432.usage", so 2*N + 7 for "/.usage" + 1 for the '\0'
    assert(filename_usage != NULL);
    sprintf(filename_usage, "%s/%s.usage", basename, basename);
    printf("Usage filename: '%s'\n", filename_usage);

    // Circuit netlist data
    uint num_inputs, num_outputs, num_gates;
    uint* inputs;
    uint* outputs;
    Gate* gates;

    // Test patterns data
    uint num_tests;
    uchar *tests_v1;
    uchar *tests_v2;
    uchar *tests_expected;

    // Fault data
    uint num_faults;
    Fault* faults;

    // Ideally, we could just now say something like:
    // read_netlist(&num_inputs, &num_outputs, &num_gates, &num_faults, &inputs, &outputs, &gates);
    // read_tests(&num_tests, &tests);
    // ------------- Begin "please move to separate function/file" section --------------

    // now we need to load in our circuit netlist
    FILE* fp = fopen(filename_input, "r");
    assert(fp != NULL);

    if (fscanf(fp, "NUMINPUTS %d\n", &num_inputs) != 1)
    {
        fprintf(stderr, "Unable to parse NUMINPUTS line!\n");
        exit(1);
    }

    printf("Detected %d inputs\n", num_inputs);
    malloc_bytes += (sizeof(uint) * num_inputs);
    inputs = (uint*) malloc(sizeof(uint) * num_inputs);
    assert(inputs != NULL);
    for (uint i = 0; i < num_inputs; i++)
    {
        if(fscanf(fp, "INPUT %d\n", &inputs[i]) != 1)
        {
            fprintf(stderr, "Unable to parse INPUT line #%d!\n", i);
            exit(1);
        }
    }


    if (fscanf(fp, "NUMOUTPUTS %d\n", &num_outputs) != 1)
    {
        fprintf(stderr, "Unable to parse NUMOUTPUTS line!\n");
        exit(1);
    }

    printf("Detected %d outputs\n", num_outputs);
    malloc_bytes += (sizeof(uint) * num_outputs);
    outputs = (uint*) malloc(sizeof(uint) * num_outputs);
    assert(outputs != NULL);
    for (uint i = 0; i < num_outputs; i++)
    {
        if (fscanf(fp, "OUTPUT %d\n", &outputs[i]) != 1)
        {
            fprintf(stderr, "Unable to parse OUTPUT line #%d!\n", i);
            exit(1);
        }
    }


    if (fscanf(fp, "NUMGATES %d\n", &num_gates) != 1)
    {
        fprintf(stderr, "Unable to parse NUMGATES line!\n");
        exit(1);
    }

    // IMPORTANT - We assume that the gates in the file are already in topological order!
    pretty_bytes(buffer, sizeof(Gate) * num_gates);
    printf("Detected %d gates (%s)\n", num_gates, buffer);
    malloc_bytes += (sizeof(Gate) * num_gates);
    gates = (Gate*) malloc(sizeof(Gate) * num_gates);
    uint type, out, in1, in2;
    for (uint i = 0; i < num_gates; i++)
    {

        if (fscanf(fp, "%d %d %d %d\n", &type, &out, &in1, &in2) != 4)
        {
            fprintf(stderr, "Unable to parse gate line #%d!\n", i);
            exit(1);
        }
        // there has to be a better way to do "if i in list of outputs"
        uchar is_output = 0;
        for (uint output_id = 0; output_id < num_outputs; output_id++)
        {
            if (i == outputs[output_id])
            {
                is_output = 1;
                break;
            }
        }
        gates[i].type = type;
        gates[i].is_output = is_output;
        gates[i].in1 = in1;
        gates[i].in2 = in2;
    }
    fclose(fp);


    // Read in tests
    fp = fopen(filename_tests, "r");
    assert(fp != NULL);

    if (fscanf(fp, "NUMTESTS %d\n", &num_tests) != 1)
    {
        fprintf(stderr, "Unable to parse NUMTESTS line!\n");
        exit(1);
    }

    uint size_v1_v2 = divide_round_up(num_inputs, 8);
    uint size_expected = divide_round_up(num_outputs, 8);
    pretty_bytes(buffer, (size_v1_v2 * 2 + size_expected) * num_tests);
    printf("Detected %d tests (%s)\n", num_tests, buffer);

    // NEW PLAN - Information stored in a compacted format, eight bits per byte, no TestPair structure, just big arrays for v1, v2, and expected, since pointers cost us 8 bytes!
    malloc_bytes += (size_v1_v2 * num_tests);
    malloc_bytes += (size_v1_v2 * num_tests);
    malloc_bytes += (size_expected * num_tests);
    tests_v1 = (uchar*) malloc(size_v1_v2 * num_tests);
    tests_v2 = (uchar*) malloc(size_v1_v2 * num_tests);
    tests_expected = (uchar*) malloc(size_expected * num_tests);
    assert(tests_v1 != NULL);
    assert(tests_v2 != NULL);
    assert(tests_expected != NULL);

    // these buffers are just for reading from the file
    malloc_bytes += (num_inputs + 1);
    char* buf_v1 = (char*) malloc(num_inputs + 1);
    assert(buf_v1 != NULL);
    malloc_bytes += (num_inputs + 1);
    char* buf_v2 = (char*) malloc(num_inputs + 1);
    assert(buf_v2 != NULL);
    malloc_bytes += (num_outputs + 1);
    char* buf_expected = (char*) malloc(num_outputs + 1);
    assert(buf_expected != NULL);
    for (uint test_id = 0; test_id < num_tests; test_id++)
    {
        if (fscanf(fp, "%s %s %s\n", buf_v1, buf_v2, buf_expected) != 3)
        {
            fprintf(stderr, "Unable to parse tests line #%d!\n", test_id);
            exit(1);
        }

        // now we need to convert the values to our new special compacted binary format
        for (uint i = 0; i < num_inputs; i++)
        {
            BIT_SET_UCHAR(tests_v1 + size_v1_v2 * test_id, i, (buf_v1[i] == '0') ? 0 : 1);
            BIT_SET_UCHAR(tests_v2 + size_v1_v2 * test_id, i, (buf_v2[i] == '0') ? 0 : 1);
        }
        for (uint i = 0; i < num_outputs; i++)
        {
            BIT_SET_UCHAR(tests_expected, size_expected * test_id + i, (buf_expected[i] == '0') ? 0 : 1);
        }
    }
    free(buf_v1);
    free(buf_v2);
    free(buf_expected);
    fclose(fp);


    // Read in list of faults
    fp = fopen(filename_faults, "r");
    assert(fp != NULL);

    if (fscanf(fp, "NUM_FAULTS %d\n", &num_faults) != 1)
    {
        fprintf(stderr, "Unable to parse NUM_FAULTS line!\n");
        exit(1);
    }

    printf("Detected %d faults\n", num_faults);
    malloc_bytes += (sizeof(Fault) * num_faults);
    faults = (Fault*) malloc(sizeof(Fault) * num_faults);
    assert(faults != NULL);
    uint net, polarity;
    for (uint ix = 0; ix < num_faults; ix++)
    {
        if (fscanf(fp, "%d %d\n", &net, &polarity) != 2)
        {
            fprintf(stderr, "Unable to parse fault line #%d!\n", ix);
            exit(1);
        }

        faults[ix].net = net;
        faults[ix].polarity = polarity;
    }
    fclose(fp);


    printf("Finished parsing files!\n");
    printf("--------------------------------------------\n");

    // ------------- End "please move to separate function/file" section --------------

    gettimeofday(&tvDoneParsing, NULL);

    /***************************************************************************
     * The Grand Plan
     * A single-fault-multiple-pattern approach:
     * 1. First, we run a parallel fault simulation to determine the fault-free
     * circuit state for every test pair.
     * 2. Then, for each fault, we determine which test pairs activate the fault
     * 3. We can ignore test pairs that do not activate the fault.
     * 4. We do another parallel fault simulation on only the patterns that
     *    activate the fault. We also have to change each state to put the X
     *    value at the necessary net before fault simulation.
    ***************************************************************************/

    /*_  ________ _____  _   _ ______ _        __
    | |/ /  ____|  __ \| \ | |  ____| |      /_ |
    | ' /| |__  | |__) |  \| | |__  | |       | |
    |  < |  __| |  _  /| . ` |  __| | |       | |
    | . \| |____| | \ \| |\  | |____| |____   | |
    |_|\_\______|_|  \_\_| \_|______|______|  |_|
    */

    // We need to store the circuit state (v1 and v2 values for all nets) for all tests.
    uint num_nets = num_inputs + num_gates;
    uint num_state_values = num_nets * 2; // times 2, since we have v1 and v2 state! see explanation above
    uint state_bytes = num_state_values + 1; // weird, but it works - Plus 1 for the extra byte at the end to track detected/not-detected.
    uint all_states_size = state_bytes * num_tests; // number of bytes required to store a state for each test
    pretty_bytes(buffer, all_states_size);
    printf("Detected %d nets, requiring %u B per state, %s total\n", num_nets, state_bytes, buffer);

    // we allocate a circuit state for each test, and then we find the fault-free values in the circuit for each test
    malloc_bytes += (all_states_size);
    uchar* fault_free_states = (uchar*) malloc(all_states_size);
    assert(fault_free_states != NULL);
    // let's set the values all to X (X=01, so 01010101, so 0x55) - This uses our new compact format for storing the state information
    memset(fault_free_states, 0x55, all_states_size);
    //print_state(fault_free_states, num_state_values);
    //print_state_raw(fault_free_states, num_state_values);

    // Set the input patterns in all fault-free states
    for (uint test_id = 0; test_id < num_tests; test_id++)
    {
        uchar* this_state = fault_free_states + test_id * state_bytes;

        // Clear the detected byte for this test
        this_state[state_bytes - 1] = 0;

        // set the input values in the fault free states
        for (uint input_id = 0; input_id < num_inputs; input_id++)
        {
            this_state[inputs[input_id] * 2]     = (BIT_GET_UCHAR(tests_v1 + size_v1_v2 * test_id, input_id) == 0 ? LOGIC_0 : LOGIC_1);
            this_state[inputs[input_id] * 2 + 1] = (BIT_GET_UCHAR(tests_v2 + size_v1_v2 * test_id, input_id) == 0 ? LOGIC_0 : LOGIC_1);
        }
    }

    // Now that we have set the input patterns in all the fault-free states, we run the simulations in parallel to find the fault-free circuit states for all tests
    // We run a parallel fault simulation to determine the fault-free circuit state for every test pair

    // the gpu needs copies of the states and the gates list-of-structs
    uchar* dev_fault_free_states;
    cuda_malloc_bytes += all_states_size;
    cudaMalloc( (void**)&dev_fault_free_states, all_states_size);
    check_cuda_errors("(cudaMalloc dev_fault_free_states)");
    assert(dev_fault_free_states != NULL);
    cudaMemcpy( dev_fault_free_states, fault_free_states, all_states_size, cudaMemcpyHostToDevice );
    check_cuda_errors("(cudaMemcpy dev_fault_free_states to GPU)");

    Gate* dev_gates;
    cuda_malloc_bytes += (sizeof(Gate) * num_gates);
    cudaMalloc( (void**)&dev_gates, sizeof(Gate) * num_gates);
    check_cuda_errors("(cudaMalloc dev_gates)");
    assert(dev_gates != NULL);
    cudaMemcpy( dev_gates, gates, sizeof(Gate) * num_gates, cudaMemcpyHostToDevice );
    check_cuda_errors("(cudaMemcpy dev_gates to GPU)");

    // Launch Kernel 1! We have a blocks with FAULTS_PER_BLOCK_KERNEL_1 threads, one thread for each test
    gettimeofday(&tvPreK1, NULL);
    cuda_fault_free_fault_sim<<< divide_round_up(num_tests, FAULTS_PER_BLOCK_KERNEL_1), FAULTS_PER_BLOCK_KERNEL_1 >>>(dev_gates, num_gates, dev_fault_free_states, state_bytes, num_tests);
    check_cuda_errors("1 (fault free fault simulation)");
    gettimeofday(&tvPostK1, NULL);
    printf("finished with fault-free responses kernel #1\n");

    cudaMemcpy( fault_free_states, dev_fault_free_states, all_states_size, cudaMemcpyDeviceToHost );
    check_cuda_errors("(cudaMemcpy dev_fault_free_states to CPU)");
#if 0
    // For TRAX Multi-Fault Injections we want to get an idea of which gates spend their time outputting 1, which increases NBTI (PMOS are turned on to output 1)
    //print_state(fault_free_states, num_state_values);
    //print_state_raw(fault_free_states, num_state_values);

    uint *count_newly_damaged_at_test = (uint*) malloc(sizeof(uint) * num_tests);
    memset(count_newly_damaged_at_test, 0, sizeof(uint) * num_tests);

    uint *net_damaged = (uint*) malloc(sizeof(uint) * num_nets);
    memset(net_damaged, 0, sizeof(uint) * num_nets);

    printf("Analyzing fault-free states to determine usage of %d nets...\n", num_nets);
    uint damage_threshold = (0.7 * 2 * num_tests);
    fp = fopen(filename_usage, "w");
    for (uint net_id = 0; net_id < num_nets; net_id++) {
        uint high_count = 0;
        printf("\r%f", (1 + net_id) / (1.0 * num_nets));
        for (uint test_id = 0; test_id < num_tests; test_id++) {
            uchar *this_state = fault_free_states + test_id * state_bytes;
            // Need to check both v1 and v2 for each test
            // In the fault-free circuit there can only be 0 and 1 values
            if (this_state[net_id * 2] == LOGIC_1) {
                high_count++;
            }
            if (this_state[net_id * 2 + 1] == LOGIC_1) {
                high_count++;
            }
            // Check if this test has caused this net to go newly-damaged
            // Accumulating damage when the gate input is low
            uint low_count = (2 * (test_id + 1)) - high_count;
            if ((!net_damaged[net_id]) && (low_count == damage_threshold)) {
                // We've seen more than 70% low counts, it's damaged
                count_newly_damaged_at_test[test_id] += 1;
                net_damaged[net_id] = 1;
            }
        }
        fprintf(fp, "%d,%f\n", net_id, high_count / (2.0 * num_tests));
    }
    printf("\n");
    fclose(fp);
    printf("Finished saving the usage info...\n");

    fp = fopen("count_newly_damaged_at_test", "w");
    for (uint test_id = 0; test_id < num_tests; test_id++) {
        fprintf(fp, "%d,%d\n", test_id, count_newly_damaged_at_test[test_id]);
    }
    fclose(fp);

    exit(1);
#endif



/*_  ________ _____  _   _ ______ _        ___
 | |/ /  ____|  __ \| \ | |  ____| |      |__ \
 | ' /| |__  | |__) |  \| | |__  | |         ) |
 |  < |  __| |  _  /| . ` |  __| | |        / /
 | . \| |____| | \ \| |\  | |____| |____   / /_
 |_|\_\______|_|  \_\_| \_|______|______| |____|
*/

    // Now, at this point, we have the fault-free responses for all tests
    // From the grand plan: "2. Then, for each fault, we determine which test pairs activate the fault."
    // Each thread corresponds with a single fault, and determines which tests activate the fault
    // NEW PLAN: We need to pack this data tighter using BIT_SET_UCHAR
    uint size_fault_activations = divide_round_up(num_tests * num_gates * 2, 8);
    pretty_bytes(buffer, size_fault_activations);
    printf("We need %s to store %ld potential fault activation bits\n", buffer, ((unsigned long) num_tests) * num_gates * 2);
    malloc_bytes += (size_fault_activations);
    uchar* fault_activations = (uchar*) malloc(size_fault_activations);
    assert (fault_activations != NULL);
    memset(fault_activations, 0, size_fault_activations);

    uchar* dev_fault_activations;
    cuda_malloc_bytes += size_fault_activations;
    cudaMalloc( (void**)&dev_fault_activations, size_fault_activations);
    check_cuda_errors("cudaMalloc (dev_fault_activations)");
    assert(dev_fault_activations != NULL);
    // The memcpy below just copies in zeros. Is there a way to get around this, maybe an initializing cudaMalloc() ? TODO
    cudaMemcpy( dev_fault_activations, fault_activations, size_fault_activations, cudaMemcpyHostToDevice );
    check_cuda_errors("cudaMemcpy (dev_fault_activations zeros to GPU)");

    // Each thread checks all tests to see which tests activate its fault.
    // dev_fault_free_states is still in the GPU, no need to copy it back and forth between kernels!
    gettimeofday(&tvPreK2, NULL);
    uint num_blocks_kernel_2 = divide_round_up(num_gates * 2, FAULTS_PER_BLOCK_KERNEL_2);
    cuda_check_fault_activations<<< num_blocks_kernel_2, FAULTS_PER_BLOCK_KERNEL_2 >>>(dev_gates, dev_fault_free_states, state_bytes, num_tests, dev_fault_activations, num_gates * 2);
    check_cuda_errors("2 (fault activations)");
    gettimeofday(&tvPostK2, NULL);

    cudaMemcpy( fault_activations, dev_fault_activations, size_fault_activations, cudaMemcpyDeviceToHost );
    check_cuda_errors("post-2 (cudaMemcpy fault_activations to CPU)");


/*_  ________ _____  _   _ ______ _        ____
 | |/ /  ____|  __ \| \ | |  ____| |      |___ \
 | ' /| |__  | |__) |  \| | |__  | |        __) |
 |  < |  __| |  _  /| . ` |  __| | |       |__ <
 | . \| |____| | \ \| |\  | |____| |____   ___) |
 |_|\_\______|_|  \_\_| \_|______|______| |____/

     * 4. We do another parallel fault simulation on only the patterns that
     *    activate the fault. We also have to change each state to put the X
     *    value at the necessary net before fault simulation (sequential?). */
    // At this point each fault has some number of tests that activate the fault.
    // These are the only tests we need to further simulate.
    // NEW PLAN: We'll definitely have > 512 activations for larger circuits, so we have two options:
    // 1. Separate kernel invocations for each fault (THIS IS THE CHOICE I MADE GOING FORWARD)
    //    + Re-use the faulty states memory for each fault (don't need to allocate gigs and gigs up front)
    //    + Can take advantage of "later faults skip most data" speedups
    //    - Kernel overhead of having thousands of kernel invocations
    // 2. One gigantic kernel
    //    - All threads must process entire circuit (so they stay in sync) and lose the "later faults skip most data" speedup
    //    - Need to allocate num_activations * state_bytes bytes of memory (which can be huge)
    //    + Only one kernel invocation, so we avoid any/all overhead with kernel calls
    // As noted above, we decided to go with option 1, which seems to be working well for now.
    // UPDATE MAY 2017 - This was possibly a bad plan. There is a lot of overhead in each kernel invocation. Going to try option 2 now.

    // Need an array of how many tests need to be run (for each fault)
    malloc_bytes += (sizeof(uint) * num_gates * 2);
    uint* num_fault_activations = (uint*) malloc(sizeof(uint) * num_gates * 2);
    assert(num_fault_activations != NULL);

    // This is the array of offsets into the big state table (for each fault)
    malloc_bytes += (sizeof(uint) * num_gates * 2);
    uint* fault_activations_offset = (uint*) malloc(sizeof(uint) * num_gates * 2);
    assert(fault_activations_offset != NULL);

    ulong total_activations = 0;
    ulong max_num_activations = 0;
    for (uint fault_id = 0; fault_id < num_gates * 2; fault_id++)
    {
        fault_activations_offset[fault_id] = total_activations;

        uint count = 0;
        for (uint test_id = 0; test_id < num_tests; test_id++)
        {
            count += BIT_GET_UCHAR(fault_activations, fault_id * num_tests + test_id);
        }
        num_fault_activations[fault_id] = count;
        //printf("Fault %8d activated by %8d tests\n", fault_id, count);

        total_activations += count;
        if (count > max_num_activations)
        {
            max_num_activations = count;
        }
    }
    printf("Max num activations: %ld\n", max_num_activations);
    //printf("------------------------------------------------\n");
    printf("total_activations: %ld\n", total_activations);

    uint size_faulty_states = state_bytes * total_activations;
    pretty_bytes(buffer, size_faulty_states);
    printf("We require %s for our faulty states!\n", buffer);
    malloc_bytes += (size_faulty_states);
    uchar* faulty_states = (uchar*) malloc(size_faulty_states);
    assert(faulty_states != NULL);

    // let's make an array of the test_id values for each fault, in order for fault_0, then fault_1, etc
    // Note, we can't merge this pair of loops with the very similar pair of loops above, because we need to know total_activations before we can malloc here
    //      It's not a big deal because this part doesn't take much of the time
    malloc_bytes += (sizeof(uint) * total_activations);
    uint* activating_test_ids = (uint*) malloc(sizeof(uint) * total_activations);
    assert(activating_test_ids != NULL);
    uint array_index = 0;
    for (uint fault_id = 0; fault_id < num_gates * 2; fault_id++)
    {
        for (uint test_id = 0; test_id < num_tests; test_id++)
        {
            if (BIT_GET_UCHAR(fault_activations, fault_id * num_tests + test_id))
            {
                activating_test_ids[array_index] = test_id;

                // Copy corresponding fault-free test row into the faulty_states array
                uchar* fault_free_state = fault_free_states + state_bytes * test_id;
                uchar* faulty_state = faulty_states + state_bytes * array_index;
                memcpy(faulty_state, fault_free_state, state_bytes);
                // Activate fault in faulty_states
                uint my_gate_id = fault_id / 2;
                // 0x80 = "activated fault" so we avoid re-evaluating the gate and losing the activated fault!
                faulty_state[my_gate_id * 2 + 1] = 0x80 |
                #if USE_TRAX
                    LOGIC_X; // Activate the fault by marking an X at the fault site in v2
                #else
                    faulty_state[my_gate_id * 2]; // Activate the fault by copying-in the v1 value (infinitely delayed transition)
                #endif

                array_index += 1;
            }
        }
    }

    // May 2017 PLAN - There's too much overhead with each kernel 3 invocation, so let's try to do one gigantic parallel simulation of all the activated faults.
    // We copy fault-free test rows from the fault-free circuit stats array based on the activating_test_ids array above, then activate the faults,
    // then simulate the circuits in parallel, then copy back the faulty circuit stats, and finally generate the dictionary data.
    uchar* dev_faulty_states;
    cuda_malloc_bytes += size_faulty_states;
    cudaMalloc( (void**)&dev_faulty_states, size_faulty_states);
    check_cuda_errors("pre-3 (cudaMalloc faulty_states)");
    assert(dev_faulty_states != NULL);

    uint dict_size = divide_round_up(num_faults * num_tests, 8);
    pretty_bytes(buffer, dict_size);
    printf("We require %s for the packed dictionary data (%d faults, %d tests)\n", buffer, num_faults, num_tests);
    malloc_bytes += dict_size;
    uchar* dict = (uchar*) malloc(dict_size);
    assert(dict != NULL);
    memset(dict, 0, dict_size);

    pretty_bytes(buffer, malloc_bytes);
    printf("malloc_bytes: %ld (%s)\n", malloc_bytes, buffer);
    pretty_bytes(buffer, cuda_malloc_bytes);
    printf("cuda_malloc_bytes: %ld (%s)\n", cuda_malloc_bytes, buffer);

    struct timeval tvStep;
    gettimeofday(&tvPreK3, NULL);
    gettimeofday(&tvStep, NULL);
#define USE_CUDA_PROFILER (0)
#if USE_CUDA_PROFILER
    cudaProfilerStart();
#endif
    // Copy the activated faulty states from cpu to gpu
    cudaMemcpy(dev_faulty_states, faulty_states, size_faulty_states, cudaMemcpyHostToDevice);
    check_cuda_errors("pre-3 (cudaMemcpy faulty_states -> dev_faulty_states - CPU-to-GPU)");

    // Run kernel 3
    cudaDeviceSynchronize(); //TODO what's this for?
    uint num_blocks = divide_round_up(total_activations, FAULTS_PER_BLOCK_KERNEL_3);
    cuda_faulty_fault_sim<<< num_blocks, FAULTS_PER_BLOCK_KERNEL_3 >>>(dev_gates, num_gates, dev_faulty_states, state_bytes, total_activations);
    check_cuda_errors("3 (faulty fault sim)");
    cudaDeviceSynchronize(); //TODO what's this for?

    // Copy the faulty states from gpu to cpu
    cudaMemcpy(faulty_states, dev_faulty_states, size_faulty_states, cudaMemcpyDeviceToHost);
    check_cuda_errors("post-3 (cudaMemcpy dev_faulty_states -> faulty_states - GPU-to-CPU)");

#if USE_CUDA_PROFILER
    cudaProfilerStop();
#endif
    gettimeofday(&tvPostK3, NULL);
    printf("\n");

    // Analyze the faulty states to generate the pass/fail dictionary data
    array_index = 0;
    for (uint fault_id = 0; fault_id < num_gates * 2; fault_id++)
    {
        for (uint test_id = 0; test_id < num_tests; test_id++)
        {
            if (BIT_GET_UCHAR(fault_activations, fault_id * num_tests + test_id))
            {
                uchar* faulty_state = faulty_states + state_bytes * array_index;
                if (faulty_state[state_bytes - 1]) {
                    // fail
                    BIT_SET_UCHAR(dict, fault_id * num_tests + test_id, 1);
                }

                array_index++;
            }
        }
    }

    // Now, write the dictionary data to disk.
    fp = fopen(filename_dictionary, "w");
    for (uint fault_list_index = 0; fault_list_index < num_faults; fault_list_index++)
    {
        for (uint test_id = 0; test_id < num_tests; test_id++)
        {
            fprintf(fp, "%d", BIT_GET_UCHAR(dict, fault_list_index * num_tests + test_id));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    gettimeofday(&tvEnd, NULL);

    //printf("----------------------------------\n");
    printf("Wrote dictionary to '%s', goodbye!\n", filename_dictionary);
    printf("----------------------------------\n");
    printf("Detailed timing information:\n");

    timeval_subtract(&tvDiff, &tvDoneParsing, &tvStart);
    printf("1. %ld.%06ld Parse the input files: \n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPreK1, &tvDoneParsing);
    printf("2. %ld.%06ld Get ready for kernel 1\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPostK1, &tvPreK1);
    printf("3. %ld.%06ld Kernel 1\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPreK2, &tvPostK1);
    printf("4. %ld.%06ld Get ready for kernel 2\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPostK2, &tvPreK2);
    printf("5. %ld.%06ld Kernel 2\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPreK3, &tvPostK2);
    printf("6. %ld.%06ld Get ready for kernel 3\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvPostK3, &tvPreK3);
    printf("7. %ld.%06ld Kernel 3 loop\n", tvDiff.tv_sec, tvDiff.tv_usec);

    timeval_subtract(&tvDiff, &tvEnd, &tvPostK3);
    printf("8. %ld.%06ld Writing dictionary to file\n", tvDiff.tv_sec, tvDiff.tv_usec);

    printf("Total time:\n");
    timeval_subtract(&tvDiff, &tvEnd, &tvStart);
    printf("   %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);


    // free up the allocated CUDA memory
    cudaFree(dev_fault_free_states);
    cudaFree(dev_gates);
    cudaFree(dev_fault_activations);
    cudaFree(dev_faulty_states);

    // free up the allocated CPU memory
    // TODO double-check all these CPU-side free calls again!
    free(filename_input);
    free(filename_tests);
    free(filename_dictionary);
    free(filename_faults);
    free(inputs);
    free(outputs);
    free(gates);
    free(tests_v1);
    free(tests_v2);
    free(tests_expected);
    free(faults);

    free(fault_free_states);
    free(fault_activations);
    free(num_fault_activations);
    free(fault_activations_offset);
    free(activating_test_ids);
    free(faulty_states);
    free(dict);

    // Apparently some profiling data is transfered asyncronously so we have to call this function to wait for those transfers to finish
    cudaDeviceReset();

    return 0;
}

