#!/usr/bin/env python
#
# This is our Gate "object", but it's really more like a struct with a few extra helper functions.
# Important features of Gate are:
#   1. Tracks values for test vectors v1 and v2 (self.values is an array-like thing, so you can say self.values[0] and self.values[1] and somegate.values[0], etc)
#   1. Storing "refernce" to it's neighboring gates (fanout and fanin gates)
#   2. Self-evaluation - The evaluate() member function looks at the current fanin gates' values to determine what the gate's output values should be.

import sys

# this file/module contains the enumerated lists of the input transitions that could cause an output glitch
from trax_glitch_transitions import glitchy_transitions

# set this variable to True to enable some debug printing
debug = False

# These are the gate evaluation functions. The input_values parameter is a list of values like this: ["0", "1", "X", "H", "0", "1"]
def eval_nand(input_values):
    """ Evaluation function for nand gates, returns "0", "1", "X", or "H". """
    if "0" in input_values: # if any inputs are 0, output is 1
        retval = "1"
    else: # no controlling values, if any are "X" then output is "X"
        if "X" in input_values:
            retval = "X"
        else:
            if "H" in input_values: # all "1" and "H" -> "H"
                retval = "H"
            else:
                retval = "0" # all "1" -> "0"
    if debug:
        print "Evaluating NAND(%s) = %s" % (str(input_values), retval)
    return retval

def eval_and(input_values):
    """ Evaluation function for and gates, returns "0", "1", "X", or "H". """
    if "0" in input_values: # if any inputs are 0, output is 0
        retval = "0"
    else: # no controlling values, if any are "X" then output is "X"
        if "X" in input_values:
            retval = "X"
        else:
            if "H" in input_values: # all "1" and "H" -> "H"
                retval = "H"
            else:
                retval = "1" # all "1" -> "1"
    if debug:
        print "Evaluating AND(%s) = %s" % (str(input_values), retval)
    return retval

def eval_nor(input_values):
    """ Evaluation function for nor gates, returns "0", "1", "X", or "H". """
    if "1" in input_values: # if any inputs are 1, output is 0
        retval = "0"
    else: # no controlling values, if any are "X" then output is "X"
        if "X" in input_values:
            retval = "X"
        else:
            if "H" in input_values: # all "0" or "H" -> "H"
                retval = "H"
            else:
                retval = "1" # all "0" -> "1"
    if debug:
        print "Evaluating NOR(%s) = %s" % (str(input_values), retval)
    return retval

def eval_or(input_values):
    """ Evaluation function for or gates, returns "0", "1", "X", or "H". """
    if "1" in input_values: # if any inputs are 1, output is 1
        retval = "1"
    else: # no controlling values, if any are "X" then output is "X"
        if "X" in input_values:
            retval = "X"
        else:
            if "H" in input_values: # all "0" or "H" -> "H"
                retval = "H"
            else:
                retval = "0" # all "0" -> "0"
    if debug:
        print "Evaluating OR(%s) = %s" % (str(input_values), retval)
    return retval

def eval_xnor(input_values):
    """ Evaluation function for xnor gates, returns "0", "1", "X", or "H". """
    if "X" in input_values: # if any inputs are X, output is X
        retval = "X"
    else: # no X values
        if "H" in input_values: # any hazard values propagate a hazard
            retval = "H"
        else: # just 0 and 1 values, do the parity function
            num_1 = len([x for x in input_values if x == "1"])
            retval = "0" if (num_1 % 2 == 1) else "1"
    if debug:
        print "Evaluating XNOR(%s) = %s" % (str(input_values), retval)
    return retval

def eval_xor(input_values):
    """ Evaluation function for xor gates, returns "0", "1", "X", or "H". """
    assert len(input_values) == 2
    if "X" in input_values: # if any inputs are X, output is X
        retval = "X"
    else: # no X values
        if "H" in input_values: # any hazard values propagate a hazard
            retval = "H"
        else: # just 0 and 1 values, do the parity function
            num_1 = len([x for x in input_values if x == "1"])
            retval = "1" if (num_1 % 2 == 1) else "0"
    if debug:
        print "Evaluating XOR(%s) = %s" % (str(input_values), retval)
    return retval

def eval_not(input_values):
    """ Evaluation function for inverter, returns "0", "1", "X", or "H". """
    assert len(input_values) == 2
    assert input_values[0] == input_values[1]
    if "X" in input_values: # if any inputs are X, output is X
        retval = "X"
    else: # no X values
        if "H" in input_values: # any hazard values propagate a hazard
            retval = "H"
        else: # just 0 or 1, invert
            retval = "1" if input_values[0] == "0" else "0"
    if debug:
        print "Evaluating NOT(%s) = %s" % (str(input_values), retval)
    return retval

def eval_buf(input_values):
    """ Evaluation function for buffer, returns "0", "1", "X", or "H". """
    assert len(input_values) == 2
    assert input_values[0] == input_values[1]
    if "X" in input_values: # if any inputs are X, output is X
        retval = "X"
    else: # no X values
        if "H" in input_values: # any hazard values propagate a hazard
            retval = "H"
        else: # just 0 or 1
            retval = input_values[0]
    if debug:
        print "Evaluating BUF(%s) = %s" % (str(input_values), retval)
    return retval

eval_functions = {}
eval_functions["0"] = eval_and
eval_functions["1"] = eval_nand
eval_functions["2"] = eval_or
eval_functions["3"] = eval_nor
eval_functions["6"] = eval_xor
eval_functions["7"] = eval_xnor
eval_functions["4"] = eval_buf
eval_functions["5"] = eval_not

# The Gate class - Mostly just a structure, but the self printing functions are useful, plus the evalute() function contains a lot of application logic. It might be a good idea to remove the extract function into the other file...
class Gate:

    # in python, the self parameter is similar to the c/c++/java "this"
    def __init__(self, gtype, outnet):
        self.gtype = gtype
        self.outnet = outnet
        self.fanin = []
        self.fanout = []
        self.values = ["X", "X"] # these values are for test vector v1 and v2
        self.fault_site = False # boolean: Is there a fault here?
        self.fault_rising = False # boolean: Is the fault of type slow-to-rise?
        self.innets = []

    def __repr__(self):
        return self.__str__()
    def __str__(self):
        retval = "Gate(gtype: '%s', outnet: '%s')\n" % (self.gtype, self.outnet)
        retval += "  Fanin gates: " + ", ".join(self.fanin) + "\n"
        retval += "  Fanout gates: " + ", ".join(self.fanout) + "\n"
        retval += "  Values: " + str(self.values) + "\n"
        return retval
    # this is a useful terse string for debugging the state of the circuit
    def short_str(self):
        return "Gate(type: %s, outnet: %s, innets: %s, values: %s)" % (self.gtype, self.outnet, str(self.innets), self.values)

    def evaluate(self, state):
        """ Look up my fanin values to see if my value can be determined. Returns a tuple of the (v1, v2) values, but leaves it up to the caller to actually update v1 and v2, so the caller can detect when things change and alert my fanout gates. """
        input_values_v1 = []
        input_values_v2 = []
        for g in self.fanin:
            input_values_v1.append(state[g].values[0])
            input_values_v2.append(state[g].values[1])
        v1 = eval_functions[self.gtype](input_values_v1)
        v2 = eval_functions[self.gtype](input_values_v2)
        # these two variables v1 and v2 are our expected output values from this gate

        # now we need to detect if v2 should result in a hazard due to this gate's inputs
        iv1 = "".join(input_values_v1)
        iv2 = "".join(input_values_v2)
        #trax_gtype = self.gtype.lower() + str(len(input_values_v1))
        if self.gtype == "0":
            trax_gtype = "and2"
        if self.gtype == "1":
            trax_gtype = "nand2"
        if self.gtype == "2":
            trax_gtype = "or2"
        if self.gtype == "3":
            trax_gtype = "nor2"
        if self.gtype == "4":
            trax_gtype = "buf1"
        if self.gtype == "5":
            trax_gtype = "inv1"
        if self.gtype == "6":
            trax_gtype = "xor2"
        if self.gtype == "7":
            trax_gtype = "xnor2"

        if debug:
            print "Input hazard detection, trax gtype: %s, v1: %s, v2: %s" % (trax_gtype, iv1, iv2)
        if trax_gtype not in glitchy_transitions: # hmm, a gate we don't know how to handle glitches with...
            if trax_gtype not in ["not1", "not", "inv1", "inv", "buf1", "buf"]: # these are gates we know we can skip, so if our mystery gate is unknown but not one of these skippable ones... we have a problem!
                print "Input hazard detection, trax gtype: %s, v1: %s, v2: %s" % (trax_gtype, iv1, iv2)
                print "Gate type '%s' not in our dictionary of glitchy transitions! Please evaluate if this is a problem." % trax_gtype
                sys.exit(1)
        else:
            for gt in glitchy_transitions[trax_gtype]:
                # iv1 and iv2 are "".join(input values), so things like "01", "1H", etc
                if gt[0] == iv1 and gt[1] == iv2:
                    if debug:
                        print "  glitchy transition match:", gt
                    v2 = "H"

        # now we need to detect if we are a fault site, and if we should activate
        if self.fault_site:
            # we can activate in three ways:
            # 1. normal gate output transition
            # 2. direct hazard activation (inputs could cause output glitch)
            # 3. indirect hazard activation (upstream gate feeds us a hazard that propagates through us)
            # Fortunately, #2 and #3 both would result in an H as our V2 output (already determined above) in the absence of activation, so we can just look for that and convert the "H" to "X"!
            if ( (    self.fault_rising and v1 == "0" and v2 == "1") or # case 1, slow-to-rise fault
                 (not self.fault_rising and v1 == "1" and v2 == "0") or # case 1, slow-to-fall fault
                 (v2 == "H") ):                                         # cases 2 and 3
                if debug:
                    print "  Fault activated, producing X value! (v1: %s, v2: %s)" % (v1, v2)
                v2 = "X"

        return (v1, v2)

