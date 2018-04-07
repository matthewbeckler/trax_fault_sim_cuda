#!/usr/bin/env python
#
# A reference implementation for the 0,1,X,H fault simulation

import sys, os, copy, time
import MLButil

# this imports the Gate class and code from the Gate.py file
from Gate import Gate

def state_dump(state):
    print "---------------------------------------------------"
    for outnet in sorted(state.keys()):
        print state[outnet].short_str()
    print "---------------------------------------------------"

def simulate(state, v1, v2, external_events):
    """ Do the actual fault simulation. This is the core of the code. Updates the circuit state in place. """

    # This is our queue of gates that need to check themselves, due to one of their inputs changing in either v1 or v2
    # We just store the output net and use the state to look it up
    events = set()

    # add input values to the queue
    test_vectors = (v1, v2)
    for i in (0, 1):
        for outnet, new_value in test_vectors[i]:
            g = state[outnet]
            if g.values[i] == new_value:
                # this input is not changing, just ignore it
                continue
            # otherwise update our value...
            g.values[i] = new_value
            # ... and tell all our fanout receivers about this change, so they can re-evaluate their value
            for fo in g.fanout:
                events.add(fo)

    # if any events were provided externally, let's add them now
    for e in external_events:
        events.add(e)

    # process the event queue
    while len(events) > 0:
        e = events.pop()
        #print "--- Updating gate with outnet =", e
        g = state[e]
        new_values = g.evaluate(state)
        if new_values[0] != g.values[0] or new_values[1] != g.values[1]:
            #print "This gate output has changed in either v1 or v2 or both!"
            g.values[0] = new_values[0]
            g.values[1] = new_values[1]
            # now tell our fanout receivers that they need to update
            for fo in g.fanout:
                events.add(fo)


import gc
def dump_garbage():
    """
    show us what's the garbage about
    """
        
    # force collection
    print "\nGARBAGE:"
    gc.collect()

    print "\nGARBAGE OBJECTS:"
    for x in gc.garbage:
        s = str(x)
        if len(s) > 80: s = s[:80]
        print type(x),"\n  ", s


def main():

    gc.enable()
    gc.set_debug(gc.DEBUG_LEAK)
    dump_garbage()

    if len(sys.argv) < 2:
        print "Usage: %s basename" % sys.argv[0]
        sys.exit(1)
    basename = sys.argv[1]

    os.chdir(basename)

    inputs = []
    outputs = []
    gates = []
    for line in file(basename + ".easy"): # easy refers to this new easy-to-parse file format, instead of nasty verilog needing a complex parser
        if line.startswith("#"):
            continue

        # the python code doesn't need the NUM* lines, but they should make it easier for the c implementation
        if line.startswith("NUMINPUTS"):
            num_inputs = int(line.strip().split()[1])
        elif line.startswith("NUMOUTPUTS"):
            num_outputs = int(line.strip().split()[1])
        elif line.startswith("NUMGATES"):
            num_gates = int(line.strip().split()[1])
        elif line.startswith("INPUT"):
            input, net = line.strip().split()
            inputs.append(net)
        elif line.startswith("OUTPUT"):
            output, net = line.strip().split()
            outputs.append(net)
        else:
            items = line.strip().split()
            gtype = items[0]
            outnet = items[1]
            innets = items[2:]
            #print "Parsed gate:", (gtype, outnet, innets)
            gates.append( (gtype, outnet, innets) )

    print "Loaded %d inputs, %d outputs, and %d gates from '%s.easy'" % (len(inputs), len(outputs), len(gates), basename)

    # we store the circuit structure and state in a hashtable that maps output net -> gate
    # note: we treat circuit inputs as specialized gates of type INPUT
    gates_by_outnet = {}

    # let's convert the gates into our structured netlist graph
    # First, we have to create all the Gate objects for the inputs and gates...
    for outnet in inputs:
        gtype = "INPUT"
        g = Gate(gtype, outnet)
        gates_by_outnet[outnet] = g
    for gtype, outnet, innets in gates:
        # gates look like this: ("NAND", "OUTNET", ["IN1", "IN2", ...])
        g = Gate(gtype, outnet)
        g.innets = innets
        gates_by_outnet[outnet] = g

    # Now that we have created all our gates, we have to go back and link them together
    for outnet in gates_by_outnet.keys():
        # for each gate, we need to set up the fanin and fanout lists
        # each gate source only knows its fanin (ie. Y = AND(A, B) means we know our fanin gates are A and B), so we set that, and then update the fanin gates' fanout lists
        g = gates_by_outnet[outnet]
        if g.gtype == "INPUT":
            # move along, nothing to do here (inputs have no fanin, our fanout list will be updated by the gates we fanout to)
            continue
        for innet in g.innets: # g.innets is a list of the nets that are input to this gate
            input_gate = gates_by_outnet[innet]
            g.fanin.append(innet) # add this gate's outnet to our fanin list
            input_gate.fanout.append(g.outnet) # add my outnet to the input's fanout list

    # read in the test vectors: (v1, v2, resp2)
    test_vectors = []
    for line in file(basename + ".tests"):
        if line.startswith("NUMTESTS"):
            continue # we don't need to parse this line
        test_vectors.append(line.strip().split())
    print "Read %d tests from '%s.tests'" % (len(test_vectors), basename)

    # load the fault list from basename.faults
    faults = [] # (outnet, rising=True or False)
    for line in file(basename + ".faults.gpu"):
        if line.startswith("NUM"):
            continue # we don't need to parse this line (only for GPU code)
        net, pol = line.strip().split()
        polarity = True if pol == "1" else False
        #print (net, polarity)
        faults.append( (net, polarity) )

    dump_garbage()
    # we need to simulate each test vector pair in turn
    responses_by_test = [] # this is where we store the matrix of responses, first index is test, second index is fault
    print "Starting fault simulation for each test, for each fault"
    time_start_ms = int(round(time.time() * 1000))
    test_id = 0
    time_start = time.time() # seconds since the epoch
    time_step = time_start
    num_tests = len(test_vectors)
    state = None
    for test in test_vectors:
        # fault simulate the circuit to get the baseline responses to v1 and v2
        v1 = zip(inputs, test[0])
        v2 = zip(inputs, test[1])
        if state:
            del state
        state = copy.deepcopy(gates_by_outnet) # WTF TODO why do we have to start over with a fresh circuit state!?!?
        simulate(state, v1=v1, v2=v2, external_events=[])
        #state_dump(state)

        #print "  Analyzing %d faults..." % len(faults)
        responses = []
        prior_outnet = None
        fault_id = 0
        for outnet, rising in faults:
            #print "Fault %d" % fault_id
            #print "-----------------------------------------------------"
            #print "    Outnet '%s', rising fault? %s, v1: %s, v2: %s" % (outnet, "true" if rising else "false", state[outnet].values[0], state[outnet].values[1])

            # find our faulty gate and mark it as the fault site of the proper polarity
            fs = state[outnet]
            fs.fault_site = True
            fs.fault_rising = rising

            if prior_outnet:
                external_events = [prior_outnet, outnet]
            else:
                external_events = [outnet]
            simulate(state, v1=[], v2=[], external_events=external_events)
            # record output
            output_v1 = "".join([state[x].values[0] for x in outputs])
            output_v2 = "".join([state[x].values[1] for x in outputs])
            resp = ""
            for i in range(len(outputs)):
                val = state[outputs[i]].values[1]
                if val != "H":
                    resp += val
                else: # we replace hazard values with the original value expected
                    resp += test[2][i]
            responses.append(resp)

            # if we aren't doing the expensive deepcopy each time, let's try this
            fs.fault_site = False
            prior_outnet = outnet

            fault_id += 1

        responses_by_test.append(responses)

        time_so_far = time.time() - time_start # seconds
        progress = (test_id + 1) / float(num_tests) # we add one here because we just finished the current test_id
        estimated_time_left = MLButil.estimated_time_left(time_so_far, progress)
        print "test_id: %5d (%3.2f%%\t - %s step, %s total, %s left)" % (test_id, 100 * progress, MLButil.pretty_time(time.time() - time_step), MLButil.pretty_time(time_so_far), MLButil.pretty_time(estimated_time_left))
        time_step = time.time()

        test_id += 1
        dump_garbage()

    time_end_ms = int(round(time.time() * 1000))
    print "----------------------------------------------------"
    print "Execution time of fault simulation code: %f seconds" % ((time_end_ms - time_start_ms) / 1000.0)
    #print "%d Fault sites:" % len(faults), " ".join(["%s/%s" % (x[0], "STR" if x[1] else "STF") for x in faults])
    #for i in range(len(responses_by_test)):
        #print "Test %4d, responses:" % i, responses_by_test[i]

    filename = basename + ".dictionary.full"
    print "Writing full response dictionary to '%s'" % filename
    with open(filename, "w") as fid:
        for fault_id in range(len(faults)):
            fid.write(" ".join([responses_by_test[test_id][fault_id] for test_id in range(len(responses_by_test))]) + " \n")
            #print "Fault %d: " % fault_id + " ".join([responses_by_test[test_id][fault_id] for test_id in range(len(responses_by_test))])

    filename = basename + ".dictionary.pf"
    print "Writing pass/fail dictionary to '%s'" % filename
    with open(filename, "w") as fid:
        for fault_id in range(len(faults)):
            fid.write("%s\n" % "".join([("1" if "X" in responses_by_test[test_id][fault_id].upper() else "0") for test_id in range(len(responses_by_test))]))


if __name__ == "__main__":
    main()

