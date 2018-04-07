#!/usr/bin/env python
#
# Prepares a circuit to run in the GPU fault simulator.
#
# Important things this script must do:
# ------------------------------------------
# 1. Gates must be topologically ordered!
# 2. Nets must be numbered starting with the output of gate 0, ending with the output of gate N-1!
# 3. Input nets will start with number N
#
# Input files:
#   basename.v
#   basename.faults
#   basename.tests (not used here, but used by the fault simulators)
# Output files:
#   basename.easy
#   basename.faults.gpu
#   basename.net_mapping

import VerilogParser
import sys, os, operator

def merge(net_pair):
    if net_pair[1] != '':
        return net_pair[0] + "_" + net_pair[1]
    else:
        return net_pair[0]


# We need to map a more loosely-defined gate name into a standardized capitalization and gate typing (NOT instead of INV, for example)
def gate_map(gate_type):
    types = ["not", "inv", "buf", "nand", "nor", "and", "or", "xnor", "xor"]
    values = ["NOT", "NOT", "BUF", "NAND", "NOR", "AND", "OR", "XNOR", "XOR"]
    gate_type_map = {"BUF": 4, "NOT": 5, "AND": 0, "NAND": 1, "OR": 2, "NOR": 3, "XOR": 6, "XNOR": 7}
    assert len(types) == len(values)
    for i in range(len(types)):
        if gate_type.startswith(types[i]):
            return gate_type_map[values[i]]
    print "Missing mapping for gate type '%s'" % gate_type
    sys.exit(1)



if len(sys.argv) < 2:
    print "Usage: %s basename" % sys.argv[0]
    sys.exit(1)
basename = sys.argv[1]
os.chdir(basename)

modules = VerilogParser.parse(basename + ".v")
assert len(modules) == 1, "This script only supports files with a single module"
mod = modules[0]

gates_by_outnet = {}
for gate in mod['gate_named']:
    gtype = gate_map(gate[0])
    outnet = merge(gate[2]['Y'])
    num_inputs = len(gate[2]) - 1
    #print gate
    assert(num_inputs <= 2), "Found a gate with more than 2 inputs! %s" % str(gate)
    in1 = merge(gate[2]['A'])
    in2 = merge(gate[2]['B']) if ('B' in gate[2]) else None
    gates_by_outnet[outnet] = [ gtype, outnet, in1, in2, None ] # None is the level
    #print "Gate, type: %d, outnet: %s, in1: %s, in2: %s" % (gtype, outnet, in1, in2)

# Let's find the level of the gates in a slow and dumb way,
# just keep iterating until we don't update any level values
done = False
while not done:
    done = True
    for outnet in gates_by_outnet:
        g = gates_by_outnet[outnet]

        if g[2] not in gates_by_outnet: # then in1 is an input
            in1_lvl = 0
        else:
            in1_lvl = gates_by_outnet[g[2]][4]
        if g[3] not in gates_by_outnet: # then in2 is an input, this also handles inv/buf situations properly
            in2_lvl = 0
        else:
            in2_lvl = gates_by_outnet[g[3]][4]

        if in1_lvl is None or in2_lvl is None:
            # we don't know enough yet, can't do anything, wait for the next pass
            done = False
        else:
            my_lvl = max(in1_lvl, in2_lvl) + 1
            g[4] = my_lvl
#for outnet in gates_by_outnet:
#    print gates_by_outnet[outnet]


# Now that we have found the level of our gates, we can easily put them in topological order
temp_gate_list = gates_by_outnet.values() # get just the gates
#print temp_gate_list
gates_in_topo_order = sorted(temp_gate_list, key=operator.itemgetter(4,3))
#print gates_in_topo_order

# we need to convert our net names into integers, starting with the gate output nets (already in topological order)
nets = []
# first the topo-sorted gate output nets
for gtype, outnet, in1, in2, level in gates_in_topo_order:
    nets.append(outnet)
# then the inputs, in the given order
for inp in mod['inputs']:
    nets.append(merge(inp))

print "Found %d nets in this circuit" % len(nets)
#print nets
net_to_index_map = {}
with open(basename + ".net_mapping", "w") as fid:
    for i in range(len(nets)):
        net_to_index_map[nets[i]] = i
        fid.write("%s\n" % nets[i])
#print net_to_index_map

# now let's read in the net names of the original gate outputs so we can convert them to net IDs
faults = []
for line in file(basename + ".faults"):
    net, polarity = line.strip().split()
    faults.append( (net_to_index_map[net], polarity) )
with open(basename + ".faults.gpu", "w") as fid:
    fid.write("NUM_FAULTS %s\n" % str(len(faults)))
    for netid, polarity in faults: # important to maintain the same ordering
        pol = 1 if polarity == "STR" else 0
        fid.write("%d %d\n" % (netid, pol) )

fid = open(basename + ".easy", "w")
fid.write("NUMINPUTS %d\n" % len(mod['inputs']))
for inp in mod['inputs']:
    fid.write("INPUT %s\n" % net_to_index_map[merge(inp)])
fid.write("NUMOUTPUTS %d\n" % len(mod['outputs']))
for outp in mod['outputs']:
    fid.write("OUTPUT %s\n" % net_to_index_map[merge(outp)])
fid.write("NUMGATES %d\n" % len(mod['gate_named']))
for gtype, outnet, in1, in2, level in gates_in_topo_order:
    inputs = [net_to_index_map[merge(gate[2][x])] for x in gate[2] if x != 'Y']
    fid.write("%d %d %d %d\n" % (gtype, net_to_index_map[outnet], net_to_index_map[in1], net_to_index_map[in1 if (in2 is None) else in2]))
    # in the case of buf/inv, we can just put in1 in there twice, the second copy won't actually be used to evalute anything

fid.close()

