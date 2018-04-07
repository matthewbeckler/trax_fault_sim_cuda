#!/usr/bin/env python
#
# Verilog Parser - Reads and parses a verilog file. Stores the design info
# in a dictionary-per-module, with a few fields within, such as for inputs,
# outputs, wires, assigns, and both named- and positionally-instantiated 
# gates/modules. Please let me know if you run across a verilog file that
# doesn't parse properly!
#
# TODO list:
#   Standardize handling of bus vs not-bus, perhaps default to buses, and offer a conversion function to switch to not-bus?
#
# Usage:
# import VerilogParser
# my_modules_list = VerilogParser.parse("file.v")
#
# Look at the end of this file for more usage examples
#
# Last updated January 8, 2013
# Matthew Beckler - mbeckler at cmu dot edu
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys, re, operator, copy

primitive_gates = ["and", "nand", "or", "nor", "xor", "xnor", "xorp", "xnorp", "not", "inv", "buf", "buff"]
primitive_gates_ports = ["Y", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V"]

def merge(pair):
    """ Used to merge a net-thing into a single string. """
    if pair[1] != '':
        return pair[0] + "_" + pair[1]
    else:
        return pair[0]

def _parse_iow(line):
    """ Parses an input, output, or wire line. Returns a list of (net_name, size), where size is like "[7:0]" or "[0:7]" (we don't care about the order) for a bus, or "" for a single net."""
    if re.search("^\w+\s*\[", line):
        # bus definition
        raw_bits, raw_names = re.search("\w+\s*\[([^\[]*)\]\s*(.*)\s*$", line).groups()
        names = map(str.strip, raw_names.strip().split(",")) # the names of the buses involved in this definition
        retval = []
        for name in names:
            retval.append( (name, "[%s]" % raw_bits) ) # This used to be " [%s]", why was it that way?
        return retval
    else:
        # single definition
        raw_names = re.search("\w+\s*(.*)\s*$", line).group(1)
        names = map(str.strip, raw_names.strip().split(","))
        retval = [(n, "") for n in names]
        return retval

def _parse_module(line):
    """ Parses a module definition line, returns the module name and a list of the module ports. """
    mod_name, mod_ports = re.search("module\s*([^\s\(]*)\s*\(\s*(.*)\s*\).*", line).groups()
    mod_ports.strip() # if there are trailing spaces inside the parentheses, this will help clean them out
    mod_ports = map(str.strip, mod_ports.split(","))
    return mod_name, mod_ports

def _parse_assign(line):
    """ Parses an "assign" line. Returns a list of (lhs, rhs) tuples. """
    #print "--------------------------------------------------------------------------------"
    #print line
    assigns = []
    assign, remainder = re.search("(assign)\s+(.*)", line).groups()
    while len(remainder) > 0:
        lhs, rest = re.search("([^=]*?)\s*=\s*(.*)", remainder).groups()

        if ":" in lhs: # bus definition
            lhs = re.search("(.*)\s*\[(.*)\]", lhs).groups()
        else:
            lhs = (lhs, "")

        if rest.startswith("{"):
            # using concatenation
            rhs, remainder = re.search("(\{\s*[^\}]*\s*\}),?\s*(.*)", rest).groups()
            rhs = (rhs, "")
        else:
            # just a regular net
            rhs, remainder = re.search("([^,]*),?\s*(.*)", rest).groups()

            if ":" in rhs: # bus definition
                rhs = re.search("(.*)\s*\[(.*)\]", rhs).groups()
            else:
                rhs = (rhs, "")

        #print "   '%s' = '%s'" % (lhs, rhs)
        assigns.append( (lhs, rhs) )
    return assigns

def _parse_instantiation(line):
    """ Parses an instantation line (module or gate). Returns a list of (type, instantiation name, io list). The IO list should be a dictionary mapping from either 0, 1, 2... or 'Y', 'A', 'B' to nets of the form (base, bus) like ('asdf', '[3:0]'). """
    # There can be multiple instantiations per line, unfortunately:
    # nand2 Xo2( .A(NotA), .B(B), .Y(line2) ), Xo3( .A(NotB), .B(A), .Y(line3) ), Xo4( .A(line2), .B(line3), .Y(Y) )
    #print line
    gtype, remainder = re.search("([^\s]+)\s+(.*)", line).groups()
    #print "--------------------------------------------------------------------------------"
    #print "gtype: '%s'" % gtype
    #print "remainder: '%s'" % remainder
    instantiations = []
    while len(remainder) > 0:
        iolist = {}
        if re.search(".*\([^\)]*\(", remainder):
            #print "  using named assignment"
            iname, raw_iolist, remainder = re.search("([^\s\(]+)\s*\(\s*(.*?\)\s*?)\s*\)\s*,?\s*(.*)", remainder).groups()
            while len(raw_iolist) > 0:
                # we actually don't have to deal with concatenations here, since .A({foo[1], foo[0]}) always catches {foo[1], foo[0]} correctly (as far as I have seen)
                name, net, raw_iolist = re.search("\.([^\)]*)\s*\(\s*([^\)]*)\s*\)\s*,?\s*(.*)", raw_iolist).groups()
                iolist[name] = parse_net_thing(net)
        else:
            #print "  using positional assignment"
            iname, raw_iolist, remainder = re.search("([^\s\(]+)\s*\(\s*(.*?)\s*\)\s*,?\s*(.*)", remainder).groups()
            #print "iname: '%s'" % iname
            #print "raw_iolist: '%s'" % raw_iolist
            #print "remainder: '%s'" % remainder
            pos = 0
            while len(raw_iolist) > 0:
                if raw_iolist.startswith("{"):
                    # this one is a concatenation
                    #print "raw_iolist: '%s'" % raw_iolist
                    net, raw_iolist = re.search("\s*(\{[^\}]+\})\s*,?\s*(.*)", raw_iolist).groups()
                    #print "net: '%s'" % net
                    #print "   : '%s'" % str(parse_net_thing(net))
                    iolist[pos] = parse_net_thing(net)
                else:
                    # just a regular net
                    #print raw_iolist
                    net, raw_iolist = re.search("([^,\s]*)\s*,?\s*(.*)", raw_iolist).groups()
                    iolist[pos] = parse_net_thing(net)
                pos += 1

        #print iolist
        instantiations.append( (gtype, iname, iolist) )

    return instantiations

def make_empty_module(name):
    this_mod = {}
    this_mod['name'] = name
    this_mod['ports'] = []
    this_mod['inputs'] = []
    this_mod['outputs'] = []
    this_mod['wires'] = []
    this_mod['assigns'] = []
    # these next two are for instantiations (gates or modules)
    # '_pos' is for unnamed instantiations, that are positional:
    #       "not1 mynot(out, in);"
    # '_named' is for .clk(foo) style named instantiations
    #       "not1 mynot(.A(in), .Y(out));"
    this_mod['gate_pos'] = []
    this_mod['gate_named'] = []
    return this_mod

def parse(filename, dump_cleaned = False):
    # read in the file
    verilog = []
    for line in file(filename):
        line = line.strip()
        if not line:
            continue
        if line.startswith("`"):
            continue
        if line.startswith("//"):
            continue
        if "//" in line:
            line = line[:line.find("//")]
        verilog.append(line)
    
    # get rid of multiple-line statements and comment lines
    verilog = " ".join([x for x in verilog if not x.startswith("//")]) # TODO do we still need this // check if we have probably handled it above?
    # get rid of multi-line comments
    verilog = re.sub("\/\*.*?\*\/", "", verilog)
    # split the text into statements
    verilog = verilog.split(";")
    verilog = map(str.strip, verilog)

    # since verilog doesn't require a semicolon after 'endmodule', it gets joined with the next line
    # so we have to split it up manually
    changed = True # since we are changing the verilog list, we want to restart the search each time we change something
    while changed:
        changed = False
        for ix in range(len(verilog)):
            if verilog[ix].startswith("endmodule"):
                # there are lots of lines that start with 'endmodule', but not all (like previously-converted ones) have a module (foo) part too
                match = re.search("module\s+(.*)", verilog[ix])
                if match:
                    nextline = match.group(1)
                    verilog[ix] = "endmodule"
                    if ix < len(verilog):
                        verilog.insert(ix + 1, nextline)
                    else:
                        verilog.append(nextline)
                    changed = True
                    break # start over in the outer while loop

    if dump_cleaned:
        # here we can dump the raw, cleaned verilog text to file, for debugging purposes:
        # sort of like looking at the output of a compiler's preprocessor
        with open(filename + ".cleaned.v", "w") as fid:
            fid.write("\n".join(verilog) + "\n")

    modules = []
    this_mod = None

    for line in verilog:
        line = line.strip()
        if line.startswith("module") and this_mod == None:
            name, ports = _parse_module(line)
            this_mod = {}
            this_mod['name'] = name
            this_mod['ports'] = ports
            this_mod['inputs'] = []
            this_mod['outputs'] = []
            this_mod['wires'] = []
            this_mod['assigns'] = []
            # these next two are for instantiations (gates or modules)
            # '_pos' is for unnamed instantiations, that are positional, eg:
            #       "not1 mynot(out, in);"
            # '_named' is for .clk(foo) style named instantiations, eg:
            #       "not1 mynot(.A(in), .Y(out));"
            this_mod['gate_pos'] = []
            this_mod['gate_named'] = []
        elif line.startswith("endmodule"):
            modules.append(this_mod)
            this_mod = None
        elif line.startswith("input"):
            list = _parse_iow(line)
            this_mod['inputs'].extend(list)
        elif line.startswith("output"):
            list = _parse_iow(line)
            this_mod['outputs'].extend(list)
        elif line.startswith("wire"):
            list = _parse_iow(line)
            this_mod['wires'].extend(list)
        elif line.startswith("assign"):
            list = _parse_assign(line)
            this_mod['assigns'].extend(list)
        elif line.startswith("parameter"):
            # TODO make it handle parameters eventually
            print "Unhandled type: 'parameter'"
        else:
            instantiations = _parse_instantiation(line)
            for gtype, instantiation_name, iolist in instantiations:
                if 0 in iolist.keys():
                    # unnamed - positional instantiation
                    this_mod['gate_pos'].append( (gtype, instantiation_name, iolist) )
                else:
                    # named instantiation
                    this_mod['gate_named'].append( (gtype, instantiation_name, iolist) )

    return modules

def is_true_primitive(gate_name):
    return gate_name.lower() in primitive_gates
def is_primitive(gate_name):
    return reduce(operator.or_, [gate_name.lower().startswith(x) for x in primitive_gates])

def convert_to_scan(modules):
    """ This function converts a sequential (DFF-based) design to a "scan-based" design by removing the DFFs and making the DFF's IO lines into module IO as required. """
    # Rules for converting a DFF into extra inputs and outputs:
    # Input to DFF is called D, output is Q (per usual)
    # 
    # 1. Add Q to inputs
    # 2. If Q is in the outputs, remove it from the outputs
    # 3. If D is not in the inputs, add D to the outputs
    # Note that we also might run into DFFs that also have a QN output that we have to handle, too

    count = 0
    for mod in modules:

        new_gate_named = [] # This will be our new list of gate_named gates, minus all the DFFs
        # I guess we are only looking at gate_named for now...
        for gate in mod['gate_named']:
            if "latch" in gate[0]:
                print "Error: We don't know how to handle latches!"
                sys.exit(1)
            elif "dff" not in gate[0]:
                new_gate_named.append(gate)
            else:
                count += 1

                d = gate[2]['d'] if ('d' in gate[2]) else gate[2]['D']
                # if D not in inputs: add D to outputs
                if d not in mod['inputs'] and d[0] != "1'b1" and d[0] != "1'b0" and d not in mod['outputs']:
                    #print "    D %s is not already an input, adding to outputs..." % str(d)
                    mod['outputs'].append(d)
                    mod['ports'].append(d[0])


                # For cleanliness sake, remove new IO lines from wires
                if d in mod['wires']:
                    mod['wires'].remove(d)


                q = gate[2]['q'] if ('q' in gate[2]) else gate[2]['Q']
                # Add Q to the inputs
                #print "    Q \"%s\" is being added to inputs..." % q
                mod['inputs'].append(q)
                mod['ports'].append(q[0])

                # If Q in outputs, remove from outputs
                if q in mod['outputs']:
                    #print "    Q \"%s\" is in outputs, removing..." % q
                    mod['outputs'].remove(q)
                    mod['ports'].remove(q[0])

                # For cleanliness sake, remove new IO lines from wires
                if q in mod['wires']:
                    mod['wires'].remove(q)


                if 'qn' in gate[2] or 'QN' in gate[2]:
                    qn = gate[2]['qn'] if ('qn' in gate[2]) else gate[2]['QN']
                    # Add QN to the inputs
                    #print "    QN \"%s\" is being added to inputs..." % qn
                    mod['inputs'].append(qn)
                    mod['ports'].append(qn[0])

                    # If QN in outputs, remove from outputs
                    if qn in mod['outputs']:
                        #print "    QN \"%s\" is in outputs, removing..." % qn
                        mod['outputs'].remove(qn)
                        mod['ports'].remove(qn[0])

                    # For cleanliness sake, remove new IO lines from wires
                    if qn in mod['wires']:
                        mod['wires'].remove(qn)


        mod['gate_named'] = new_gate_named

    return count


def convert_to_named_ports(modules):
    """ This function converts the instantiations of gates and modules from positional to named. Modifies the module list in place! This function is so important we should almost call it automatically right after parsing. """
    # For this example, we want to convert all the gates in gate_pos to be in gate_named
    # We assume that the order of inputs is Y, A, B, C, D, E,...., that is, the output is first
    # We will only do this automatic conversion for the basic gates:
    #   and, nand, or, nor, xor, xnor, inv, not, buf
    #print "--- Converting module gates to named-instantiation ---"
    num_pos = 0
    for mod in modules:
        #print "--- Converting module '%s' to named-instantiation ---" % mod['name']
        num_pos += len(mod['gate_pos'])
        while len(mod['gate_pos']) > 0:
            gate = mod['gate_pos'][0]
            if gate[0] in [x['name'] for x in modules]:
                #print "gate %s is in our list of modules..." % gate[0]
                thismod = modules[[x['name'] for x in modules].index(gate[0])]
                connections = {}
                for index, net in gate[2].iteritems():
                    connections[thismod['ports'][index]] = net
                new_gate = (gate[0], gate[1], connections)
                mod['gate_named'].append(new_gate)
                mod['gate_pos'].pop(0)
            else:
                #print "gate %s is not in our list of modules" % gate[0]
                if is_primitive(gate[0]):
                    # nice, it's a primitive
                    connections = {}
                    for index, net in gate[2].iteritems():
                        connections[primitive_gates_ports[index]] = net
                    new_gate = (gate[0], gate[1], connections)
                    mod['gate_named'].append(new_gate)
                    mod['gate_pos'].pop(0)
                elif "dff" in gate[0]:
                    # dff, a bit special - TODO does this handle the clock signal correctly?
                    new_gate = ("dff", gate[1], {'Q': gate[2][0], 'D': gate[2][1]})
                    print "new gate from dff", gate
                    print "    ", new_gate
                    mod['gate_named'].append(new_gate)
                    mod['gate_pos'].pop(0)
                else:
                    print "Weird, gate '%s' is not in our module list and is also not a primitive. Exiting." % gate[0]
                    exit()
    return num_pos

def bus_range(bus_range_string):
    """ Creates an iterable range based on the bus string, something like " [2:0]". """
    s = bus_range_string.strip()
    if "[" in s:
        s = s[1:-1]
        a, b = map(int, s.split(":"))
    else:
        a, b = map(int, s.split(":"))
    if a < b:
        return range(a, b + 1)
    else:
        return range(a, b - 1, -1)

def convert_to_no_buses_iow(pair_list):
    """ Pass in a list of (name, bus) pairs. Handles boring cases like "net5" -> "net5", single nets that are part of a bus like "net[5]" -> "net_5", and more complex cases of ("net", "[3:0]") -> [("net_3", ''), ("net_2", ''), ("net_1", ''), ("net_0", '')], and situations with concatenations like this ("{net[1], net[0]}", '') -> [("net_1", ""), ("net_0", "")]. """
    new_items = []
    for name, bus in pair_list:
        if bus == "":
            if "{" in name: # drat, a concatenation!
                rest = name.strip()[1:-1] # trim off the leading '{' and trailing '}' characters
                while rest != "":
                    net, rest = map(str.strip, re.search("\s*([^,]+)\s*,?\s*(.*)\s*", rest).groups())
                    if ":" in net: # bus range
                        base, start, end = re.search("\s*([^\s\[]+)\s*\[\s*(\d+)\s*:\s*(\d+)\s*\].*", net).groups()
                        for i in bus_range("[%s:%s]" % (start, end)):
                            new_items.append( (base + "_" + str(i), "") )
                    else:
                        new_items.append( (net_bus_conv(net), "") )
            else: # not a concatenation
                new_items.append( (net_bus_conv(name), "") )
        else:
            for i in bus_range(bus):
                new_items.append( (name + "_" + str(i), "") )
    return new_items

def regularize_assign_side(p):
    if p[1] == "":
        if "[" in p[0]: # single assignment with a single net from a bus
            return "_".join(re.search("([^\[]+)\[([^\]]+)\]", p[0]).groups())
        else: # not a bus at all
            return p[0]
    else: # whole bus assignment
        if "[" in p[0]:
            print "Weird, this assign is a bus but also has '[' in the net name, bailing...",  str(p)
            sys.exit(1)
        else: # assigning a whole bus
            new_items = []
            for i in bus_range(p[1]):
                new_items.append(p[0] + "_" + str(i))
            return new_items

def net_bus_conv(net):
    """ Takes a net name like "asdf[3]" and makes it into "asdf_3". """
    if "[" in net:
        return "_".join(re.search("\s*([^\s\[]+)\s*\[\s*([^\s\]]+)\s*\].*", net).groups())
    else:
        return net

def parse_net_thing(net):
    """ Parses things like "net5", "net[5]", and "net[5:0]" to the usual (net, width) pair. """
    if "{" in net: # concatenation, drat!
        return (net, '') # we don't touch it
    else:
        if ":" in net: # this is a bus range
            base, start, end = re.search("\s*([^\s\[]+)\s*\[\s*(\d+)\s*:\s*(\d+)\s*\].*", net).groups()
            return (base, "[%s:%s]" % (start, end))
        else: # just a plain old net like "net5" or a single bit from a bus like "net[5]"
            return (net, '')

def convert_instantiations(modules, mod_name, old_modules):
    """ Convert instantiations for the specified module, ONLY WORKS FOR gate_named GATES, updates the module in place. """

    #print "----------------------------------------------------------- A"
    #print "convert_assignments called"
    #print "    mod_name: '%s'" % mod_name

    # extract our module from the big list of modules
    my_mod = [x for x in modules if x['name'] == mod_name]
    assert len(my_mod) == 1
    my_mod = my_mod[0]
    
    if len(my_mod['gate_pos']) > 0:
        print "WARNING: Calling convert_instantiations on a module with positionally-instantiated gates and modules is not supported. Be sure to call the convert_to_named_ports function first!. We will ignore these positionally-instantiated gates, but it's not going to work anyway."

    # here's the plan: We go through and build up a list of new_gates, and then replace the old gates
    new_gates = []
    for type, name, cons in my_mod['gate_named']:
        #print "    Looking at instantiation:"
        #print "     ", (type, name, cons)
        new_cons = {}

        # this instantiation might be another module, or might be a primitive gate, let's check
        if type in [x["name"] for x in modules]: # we're instantiating another module!
            mod = [x for x in old_modules if x["name"] == type]
            assert len(mod) == 1
            mod = mod[0]
            #print "This instantiation is another module:", mod['name']
            #print "   ", mod['inputs']
            #print "   ", mod['outputs']
            #print (type, name, cons)
            for dest, src in cons.iteritems():
                #print "dest: '%s', src: '%s'" % (dest, src)
                # we need to determine the sizes of both sides of the instantiation
                net_dest = None

                temp = [x for x in mod['inputs'] if x[0] == dest]
                if len(temp) > 0:
                    net_dest = temp[0]
                    #print "'%s' is an input to the sub-module.  Inside the module it is: %s" % (dest, net_dest)

                temp = [x for x in mod['outputs'] if x[0] == dest]
                if len(temp) > 0:
                    net_dest = temp[0]
                    #print "'%s' is an output from the sub-module. Inside the module it is: %s" % (dest, net_dest)

                if net_dest is None: # I'm not sure if this can even happen, but we definitely should check it!
                    print "Warning: net_dest is still equal to None, so I guess it wasn't an input nor an output?"
                    exit()

                # now to figure out the size of the source (in this current module)
                net_src = None
                src_base = src[0]
                if "[" in src_base: # this handles the case where the source net is like "asdf[0]"
                    src_base = src_base.split("[")[0]

                temp = [x for x in my_mod['inputs'] if x[0] == src_base]
                if len(temp) > 0:
                    net_src = temp[0]
                    #print "'%s' is an input to this module. Here it's called: %s" % (src, net_src)

                temp = [x for x in my_mod['outputs'] if x[0] == src_base]
                if len(temp) > 0:
                    net_src = temp[0]
                    #print "'%s' is an output from this module. Here it's called: %s" % (src, net_src)

                temp = [x for x in my_mod['wires'] if x[0] == src_base]
                if len(temp) > 0:
                    net_src = temp[0]
                    #print "'%s' is a local declared wire in this module. Here it's called: %s" % (src, net_src)

                if net_src is None:
                    #print "Warning: net_src is still equal to None (src = '%s'), so I guess it wasn't an input nor an output nor a declared wire. We assume it is an undeclared scalar (single bit) wire." % str(src)
                    net_src = src
                    #print src

                # now we have figured out the source and destination information
                #print "source: %s, destination: %s" % (net_src, net_dest)
                if src[1] != "" or "[" in src[0]: # the bus stuff is specified here specially, so it overrides (ie, if you connect only some bits of a local bus to a sub-module)
                    net_src_expanded = convert_to_no_buses_iow([src])
                else:
                    net_src_expanded = convert_to_no_buses_iow([net_src])
                net_dest_expanded = convert_to_no_buses_iow([net_dest])
                if len(net_src_expanded) != len(net_dest_expanded):
                    print "WARNING: When converting instantiations, the expanded lists were inequal!"
                    print "module:", mod['name']
                    print "src: '%s', dest: '%s'" % (str(net_src), str(net_dest))
                    print "src: '%s', net_src: '%s'" % (src, net_src)
                    print "net_src_expanded:"
                    print "   ", net_src_expanded
                    print "net_dest_expanded:"
                    print "   ", net_dest_expanded
                    exit(1)

                for i in range(len(net_src_expanded)):
                    # for each item in both lists, we add them as new connections to new_cons
                    new_cons[net_dest_expanded[i][0]] = net_src_expanded[i] # we only take the base net for the destination name


        else:
            #print "This instantiation is a primitive gate"
            for port_name, net in cons.iteritems():
                new_cons[port_name] = (net_bus_conv("".join(net)), "")
                #print "%s : %s -> %s" % (port_name, net, new_cons[port_name])

        new_gates.append( (type, name, new_cons) )
        #print "        becomes"
        #print "     ", (type, name, new_cons)

    my_mod['gate_named'] = new_gates
    #print "----------------------------------------------------------- B"

def convert_to_no_buses(modules):
    """ This function goes through all the inputs, outputs, ports, assigns, wires, and named instantiations to "expand" the bus notation to make scan easier. """

    old_modules = copy.deepcopy(modules)
    
    for mod in modules:
        #print "----------------------------------------------------------- D"
        #print "Processing module '%s'" % mod['name']

        # deal with assignments - this is the most difficult part
        convert_instantiations(modules, mod['name'], old_modules)

        mod['inputs'] = convert_to_no_buses_iow(mod['inputs'])
        mod['outputs'] = convert_to_no_buses_iow(mod['outputs'])
        mod['wires'] = convert_to_no_buses_iow(mod['wires'])

        new_pairs = []
        for lhs, rhs in mod['assigns']:
            #print "assign lhs: '%s', rhs: '%s'" % (lhs, rhs)
            lhs_items = regularize_assign_side(lhs)
            rhs_items = regularize_assign_side(rhs)
            if type(lhs_items) == type([]) and type(rhs_items) == type([]): # two buses
                if len(lhs_items) != len(rhs_items):
                    print "Error! We had an assign with two buses that are different sizes, bailing..."
                    print "'%s' = '%s'" % (str(lhs), str(rhs))
                    sys.exit(1)
                for i in range(len(lhs_items)):
                    new_pairs.append( ((lhs_items[i], ""), (rhs_items[i], "")) )
            elif type(lhs_items) == type("asdf") and type(rhs_items) == type("asdf"): # two non-buses
                new_pairs.append( ((lhs_items, ""), (rhs_items, "")) )
            else:
                print "Error! We had an assign with one bus and one non-bus, bailing..."
                print "'%s' = '%s'" % (str(lhs), str(rhs))
                sys.exit(1)
        mod['assigns'] = new_pairs

        # regenerate the ports from the new list of inputs and outputs - TODO make this use the original ordering of ports, not just outputs then inputs (although it should be just fine after converting to all named instantiations only)
        mod['ports'] = []
        mod['ports'].extend(map(operator.itemgetter(0), mod['outputs']))
        mod['ports'].extend(map(operator.itemgetter(0), mod['inputs']))

        #print "Results:"
        #print render_verilog([mod])

def convert_to_two_fanin(modules):
    """ This function converts all the by-name-instantiated gates to have two fan-ins. Modifies the module list in place! Returns (num_fixed, num_new, list of original gate output net names). """

    num_fixed = 0
    num_new = 0
    original_gate_output_net_names = [] # this will allow the GPU code to know which faults it has to do
    new_net_names = [] # this helps us make the list above
    gate_types =  ["and", "nand", "or", "nor", "xor", "xnor"]
    inverting_gates = ["nand", "nor", "xnor"]
    module_names = [x['name'] for x in modules]
    for mod in modules:
        #print "--- Fixing fan-in for module '%s' ---" % mod['name']
        if len(mod['gate_pos']) > 0:
            print "Module '%s' has positionally-instantiated gates/modules! This function only works on named-instantiated gates. Exiting!"
            sys.exit(1)
        if len(mod['gate_named']) > 0:
            finished = False
            while not finished: # we have to do this silly wrapper loop because we'll be removing and adding gates to mod['gate_named'] and that would mess up a simple loop
                finished = True
                for gate in mod['gate_named']:
                    gtype, gname, cons = gate
                    newcons = {}
                    for key, value in cons.iteritems():
                        newcons[key] = merge(value)
                    cons = newcons
                    if gtype in module_names:
                        # it's ok, this is one of our instantiated modules, we can just skip it
                        continue
                    clean_type = re.search("^([^\d]*)\d*$", gtype).groups()[0]
                    if clean_type not in gate_types:
                        if clean_type in ["inv", "buf", "not"]:
                            # it's ok, we can just skip these
                            continue
                        else:
                            print "Found an unknown gate type '%s', exiting!" % gtype
                            sys.exit(1)
                    base_type = clean_type # we have to break nand gates into a bunch of and + one final nand gate

                    size = len(cons) - 1
                    if size <= 2: # this gate is small enough
                        if (cons['Y'], "") not in new_net_names and cons['Y'] not in original_gate_output_net_names:
                            original_gate_output_net_names.append(cons['Y'])
                        continue
                    # else, this gate needs to be split up
                    # also, it means we are not finished
                    finished = False
                    original_gate_output_net_names.append(cons['Y'])

                    if clean_type in inverting_gates:
                        invmap = {"nand": "and", "nor": "or", "xnor": "xor"}
                        pure_type = invmap[clean_type]
                    else:
                        pure_type = clean_type


                    new_gates = []
                    if size == 3:
                        # two new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['C'], ""), 'B': (cons['Y'] + "_t0", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(1)] )
                    elif size == 4:
                        # three new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t1", {'Y': (cons['Y'] + "_t1", ""), 'A': (cons['C'], ""), 'B': (cons['D'], "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['Y'] + "_t0", ""), 'B': (cons['Y'] + "_t1", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(2)] )
                    elif size == 5:
                        # four new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t1", {'Y': (cons['Y'] + "_t1", ""), 'A': (cons['C'], ""), 'B': (cons['D'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t2", {'Y': (cons['Y'] + "_t2", ""), 'A': (cons['Y'] + "_t0", ""), 'B': (cons['Y'] + "_t1", "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['E'], ""), 'B': (cons['Y'] + "_t2", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(3)] )
                    elif size == 6:
                        # five new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t1", {'Y': (cons['Y'] + "_t1", ""), 'A': (cons['C'], ""), 'B': (cons['D'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t2", {'Y': (cons['Y'] + "_t2", ""), 'A': (cons['E'], ""), 'B': (cons['F'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t3", {'Y': (cons['Y'] + "_t3", ""), 'A': (cons['Y'] + "_t0", ""), 'B': (cons['Y'] + "_t1", "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['Y'] + "_t3", ""), 'B': (cons['Y'] + "_t2", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(4)] )
                    elif size == 7:
                        # six new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t1", {'Y': (cons['Y'] + "_t1", ""), 'A': (cons['C'], ""), 'B': (cons['D'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t2", {'Y': (cons['Y'] + "_t2", ""), 'A': (cons['E'], ""), 'B': (cons['F'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t3", {'Y': (cons['Y'] + "_t3", ""), 'A': (cons['Y'] + "_t0", ""), 'B': (cons['Y'] + "_t1", "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t4", {'Y': (cons['Y'] + "_t4", ""), 'A': (cons['Y'] + "_t2", ""), 'B': (cons['G'], "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['Y'] + "_t3", ""), 'B': (cons['Y'] + "_t4", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(5)] )
                    elif size == 8:
                        # seven new gates
                        new_gates.append( (pure_type + "2", gname + "_t0", {'Y': (cons['Y'] + "_t0", ""), 'A': (cons['A'], ""), 'B': (cons['B'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t1", {'Y': (cons['Y'] + "_t1", ""), 'A': (cons['C'], ""), 'B': (cons['D'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t2", {'Y': (cons['Y'] + "_t2", ""), 'A': (cons['E'], ""), 'B': (cons['F'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t3", {'Y': (cons['Y'] + "_t3", ""), 'A': (cons['G'], ""), 'B': (cons['H'], "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t4", {'Y': (cons['Y'] + "_t4", ""), 'A': (cons['Y'] + "_t0", ""), 'B': (cons['Y'] + "_t1", "")}) )
                        new_gates.append( (pure_type + "2", gname + "_t5", {'Y': (cons['Y'] + "_t5", ""), 'A': (cons['Y'] + "_t2", ""), 'B': (cons['Y'] + "_t3", "")}) )
                        new_gates.append( (clean_type + "2", gname, {'Y': (cons['Y'], ""), 'A': (cons['Y'] + "_t4", ""), 'B': (cons['Y'] + "_t5", "")}) )
                        new_net_names.extend( [(cons['Y'] + "_t%d" % i, "") for i in range(6)] )
                    else:
                        print "Unsupported gate size? %d" % size
                        sys.exit(0)

                    mod['gate_named'].remove(gate)
                    mod['gate_named'].extend(new_gates)
                    num_fixed += 1
                    num_new += len(new_gates)
    
    return (num_fixed, num_new, original_gate_output_net_names)

def remap(mapping, original_io_dict):
    #print "remap called:"
    #print "   ", mapping
    #print "   ", original_io_dict
    retval = {}
    for io in original_io_dict:
        net = original_io_dict[io]
        if io not in mapping:
            retval[io] = net
        else:
            retval[mapping[io]] = net
    #print "   ", retval
    return retval

def lookup(dict, keys):
    """ Pass in a dictionary, and it will check each key in turn and return the value of the first key present. """
    for k in keys:
        if k in dict:
            return dict[k]
    # TODO how to handle problems?

def simplify_gates(modules):
    """ Convert complex gates into primitives ( and convert their stupid gate and port names into standard gate and port names ). Modifies the modules list IN PLACE! Returns (num_removed_complex, num_added_primitive). """

    num_complex = 0
    num_added_primitive = 0
    primitive_mapping = {"Q": "Y", "QN": "Y", "Z": "Y", "ZN": "Y", "IN1": "A", "IN2": "B", "IN3": "C", "IN4": "D", "IN5": "E", "IN6": "F", "IN7": "G", "IN8": "H", "IN9": "J", "INP": "A"}
    for mod in modules:
        assert len(mod['gate_pos']) == 0
        #print "--- Looking for complex gates in module '%s' ---" % mod['name']
        new_gate_named = []
        for gate in mod['gate_named']:
            base = gate[1]

            if gate[0] in [x['name'] for x in modules]:
                #print "this module ('%s') is defined in the file, so we will skip it" % gate[0]
                new_gate_named.append(gate)
            elif gate[0].lower().startswith("haddx"):
                # this half adder has inputs named A0 and B0 with outputs named C1 and SO
                new_gate_named.append( ("and2", base + "_c1", {"Y": gate[2]["C1"], "A": gate[2]["A0"], "B": gate[2]["B0"]}) )
                new_gate_named.append( ("xor2", base + "_s1", {"Y": gate[2]["SO"], "A": gate[2]["A0"], "B": gate[2]["B0"]}) )
                num_complex += 1
                num_added_primitive += 2
            elif gate[0].lower().startswith("faddx"):
                # this full adder has inputs named A, B, and CI with outputs named S and CO
                new_gate_named.append( ("xor3", base + "_s1", {"Y": gate[2]["S"], "A": gate[2]["A"], "B": gate[2]["B"], "C": gate[2]["CI"]}) )
                new_gate_named.append( ("and2", base + "_c1", {"Y": (base + "_c1_Y", ""), "A": gate[2]["A"], "B": gate[2]["B"]}) )
                new_gate_named.append( ("and2", base + "_c2", {"Y": (base + "_c2_Y", ""), "A": gate[2]["A"], "B": gate[2]["CI"]}) )
                new_gate_named.append( ("and2", base + "_c3", {"Y": (base + "_c3_Y", ""), "A": gate[2]["CI"], "B": gate[2]["B"]}) )
                new_gate_named.append( ("or3", base + "_c4", {"Y": gate[2]["CO"], "A": (base + "_c1_Y", ""), "B": (base + "_c2_Y", ""), "C": (base + "_c3_Y", "")}) )
                num_complex += 1
                num_added_primitive += 5
            elif gate[0].lower().startswith("ao22x"):
                # for the ao22 gate, we need to add a pair of and2 gates driving an or2 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("or2", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", "")}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("ao222x"):
                # for the ao222 gate, we need to add a trio of and2 gates driving an or3 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("and2", base + "_temp2", {"Y": (base + "_temp2_Y", ""), "A": gate[2]["IN5"], "B": gate[2]["IN6"]}) )
                new_gate_named.append( ("or3", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": (base + "_temp2_Y", "")}) )
                num_complex += 1
                num_added_primitive += 4
            elif gate[0].lower().startswith("ao221x"):
                # for the ao221 gate, we need to add a pair of and2 gates driving an or3 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("or3", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": gate[2]["IN5"]}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("ao21x"):
                # for the ao21 gate, we need to add an and2 driving an or2 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": gate[2]["IN3"]}) )
                num_complex += 1
                num_added_primitive += 2
            elif gate[0].lower().startswith("aoi21x"):
                # for the aoi21 gate, we need to add an and2 driving a nor2 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("nor2", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": gate[2]["IN3"]}) )
                num_complex += 1
                num_added_primitive += 2
            elif gate[0].lower().startswith("aoi22x"):
                # for the aoi22 gate, we need to add a pair of and2 gates driving a nor2 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("nor2", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", "")}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("aoi222x"):
                # for the aoi222 gate, we need to add a trio of and2 gates driving a nor3 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("and2", base + "_temp2", {"Y": (base + "_temp2_Y", ""), "A": gate[2]["IN5"], "B": gate[2]["IN6"]}) )
                new_gate_named.append( ("nor3", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": (base + "_temp2_Y", "")}) )
                num_complex += 1
                num_added_primitive += 4
            elif gate[0].lower().startswith("aoi221x"):
                # for the aoi221 gate, we need to add a pair of and2 gates driving a nor3 gate
                new_gate_named.append( ("and2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("nor3", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": gate[2]["IN5"]}) )
                num_complex += 1
                num_added_primitive += 3

            elif gate[0].lower().startswith("oa22x"):
                # for the oa22 gate, we need to add a pair of or2 gates driving an and2 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("and2", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", "")}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("oa221x"):
                # for the oa221 gate, we need to add a pair of or2 gates driving an and3 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("and3", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": gate[2]["IN5"]}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("oa222x"):
                # for the oa222 gate, we need to add a trio of or2 gates driving an and3 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("or2", base + "_temp2", {"Y": (base + "_temp2_Y", ""), "A": gate[2]["IN5"], "B": gate[2]["IN6"]}) )
                new_gate_named.append( ("and3", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": (base + "_temp2_Y", "")}) )
                num_complex += 1
                num_added_primitive += 4
            elif gate[0].lower().startswith("oa21x"):
                # for the oa21 gate, we need to add an or2 driving an and2
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("and2", base, {"Y": gate[2]["Q"], "A": (base + "_temp0_Y", ""), "B": gate[2]["IN3"]}) )
                num_complex += 1
                num_added_primitive += 2
            elif gate[0].lower().startswith("oai21x"):
                # for the oai21 gate, we need to add an or2 driving a nand2
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("nand2", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": gate[2]["IN3"]}) )
                num_complex += 1
                num_added_primitive += 2
            elif gate[0].lower().startswith("oai22x"):
                # for the oai22 gate, we need to add a pair of or2 gates driving a nand2 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("nand2", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", "")}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("oai221x"):
                # for the oai221 gate, we need to add a pair of or2 gates driving a nand3 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("nand3", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": gate[2]["IN5"]}) )
                num_complex += 1
                num_added_primitive += 3
            elif gate[0].lower().startswith("oai222x"):
                # for the oai222 gate, we need to add a trio of or2 gates driving a nand3 gate
                new_gate_named.append( ("or2", base + "_temp0", {"Y": (base + "_temp0_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["IN2"]}) )
                new_gate_named.append( ("or2", base + "_temp1", {"Y": (base + "_temp1_Y", ""), "A": gate[2]["IN3"], "B": gate[2]["IN4"]}) )
                new_gate_named.append( ("or2", base + "_temp2", {"Y": (base + "_temp2_Y", ""), "A": gate[2]["IN5"], "B": gate[2]["IN6"]}) )
                new_gate_named.append( ("nand3", base, {"Y": lookup(gate[2], ("Q", "QN")), "A": (base + "_temp0_Y", ""), "B": (base + "_temp1_Y", ""), "C": (base + "_temp2_Y", "")}) )
                num_complex += 1
                num_added_primitive += 4

            elif gate[0].lower().startswith("mux2"):
                # for the mux2 gate, we need to add an inverter on S, two and2 gates, and a final or2 gate
                new_gate_named.append( ("inv1", base + "_sbar",   {"Y": (base + "_sbar_Y", ""), "A": gate[2]["S"]}) )
                new_gate_named.append( ("and2", base + "_and2_0", {"Y": (base + "_and2_0_Y", ""), "A": gate[2]["IN2"], "B": (base + "_sbar_Y", "")}) )
                new_gate_named.append( ("and2", base + "_and2_1", {"Y": (base + "_and2_1_Y", ""), "A": gate[2]["IN1"], "B": gate[2]["S"]}) )
                new_gate_named.append( ("or2", base, {"Y": gate[2]["Q"], "A": (base + "_and2_0_Y", ""), "B": (base + "_and2_1_Y", "")}) )
                num_complex += 1
                num_added_primitive += 4

            elif gate[0].lower().startswith("xnor"):
                new_gate_named.append( ("xnor%d" % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("xor"):
                new_gate_named.append( ("xor%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("nor"):
                new_gate_named.append( ("nor%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("or"):
                new_gate_named.append( ("or%d"   % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("nand"):
                new_gate_named.append( ("nand%d" % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("and"):
                new_gate_named.append( ("and%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )

            elif gate[0].lower().startswith("inv"):
                new_gate_named.append( ("inv%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("not"):
                new_gate_named.append( ("inv%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )
            elif gate[0].lower().startswith("buf") or gate[0].lower().startswith("nbuf"):
                new_gate_named.append( ("buf%d"  % (len(gate[2]) - 1), gate[1], remap(primitive_mapping, gate[2])) )

            else: # the gate is not a primitive gate, and it's not one of our modules, so it's something we don't know about and need to fix!
                print "    Gate '%s' is not primitive, and not one of the defined modules! That's pretty weird, we probably need to add more code to handle this previously-unseed complex gate!" % gate[0]
                print gate
                sys.exit(1)
        mod['gate_named'] = new_gate_named

    return num_complex, num_added_primitive


def dump_modules(modules):
    print "--------------------------------------------------------------------------------"
    for mod in modules:
        dump_module(mod)
        print "--------------------------------------------------------------------------------"

def dump_module(mod):
    """ Dumps the information for the given module to stdout, not super easy to read, but it's good in a pinch. """
    print "Name:", mod['name']
    print "Ports:", mod['ports']
    print "Inputs:", mod['inputs']
    print "Outputs:", mod['outputs']
    print "Wires:", mod['wires']
    print "Assigns:", mod['assigns']
    print "Gate Pos:"
    for g in mod['gate_pos']:
        print " ", g
    print "Gate Named:"
    for g in mod['gate_named']:
        print " ", g

def write_verilog(modules, output_filename):
    """ Writes a list of modules back into a verilog file, with clean and consistent formatting. """
    with open(output_filename, "w") as fid:
        fid.write(render_verilog(modules))

def render_verilog(modules):
    """ Creates the verilog code for these modules, as a string. """
    # tmax has a line length limit of 50000 or something, but let's make nicely human readable files
    output = ""
    maxlen = 80
    for m in modules:
        start = "module %s (" % m['name']
        output += start
        linelen = len(start)
        for name in m['ports'][:-1]: # save the last one since it doesn't get a comma after it
            if linelen + len(name) >= maxlen:
                output += "\n    "
                linelen = 4
            output += name + ", "
            linelen += len(name) + 2
        # deal with the last arg here, different than the rest (no trailing comma)
        last_port = m['ports'][-1]
        if linelen + len(last_port) >= maxlen:
            output += "\n    "
        output += m['ports'][-1] + ");\n\n"

        for name, bus in m['inputs']:
            output += "    input %s %s;\n" % (bus, name)
        output += "\n"

        for name, bus in m['outputs']:
            output += "    output %s %s;\n" % (bus, name)
        output += "\n"

        if len(m['wires']) > 0:
            for name, bus in m['wires']:
                output += "    wire %s %s;\n" % (bus, name)
            output += "\n"

        if len(m['assigns']) > 0:
            for lhs, rhs in m['assigns']:
                output += "    assign %s = %s;\n" % ("".join(lhs), "".join(rhs))
            output += "\n"


        if len(m['gate_pos']) > 0:
            for gp in m['gate_pos']:
                cons = []
                for index in range(len(gp[2])):
                    cons.append("".join(merge(gp[2][index])))
                output += "    %s %s(%s);\n" % (gp[0], gp[1], ", ".join(cons))
            output += "\n"

        if len(m['gate_named']) > 0:
            for gn in m['gate_named']:
                cons = []
                # we want to put the output connections first, for increased compatability with Osei's v2i script
                if is_primitive(gn[0]):
                    if "Y" in gn[2]:
                        cons.append(".Y(%s)" % (merge(gn[2]['Y'])))
                    if "Z" in gn[2]:
                        cons.append(".Z(%s)" % (merge(gn[2]['Z'])))
                    if "ZN" in gn[2]:
                        cons.append(".ZN(%s)" % (merge(gn[2]['ZN'])))
                    for c in gn[2]:
                        if c not in ("Y", "Z", "ZN"):
                            cons.append(".%s(%s)" % (c, merge(gn[2][c])))
                else: # this is a module instantiation
                    the_mod = None
                    for mod in modules:
                        if mod['name'] == gn[0]:
                            the_mod = mod
                            break
                    for out in the_mod['outputs']:
                        cons.append(".%s(%s)" % (out[0], merge(gn[2][out[0]])))
                    for inp in the_mod['inputs']:
                        cons.append(".%s(%s)" % (inp[0], merge(gn[2][inp[0]])))
                output += "    %s %s(%s);\n" % (gn[0], gn[1], ", ".join(cons))
            output += "\n"

        output += "endmodule /* %s */\n\n" % m['name']
    return output

def main():
    import sys
    if len(sys.argv) < 2:
        print "Usage: %s verilog_filename" % sys.argv[0]
        print "Running this script will parse a verilog file, report any errors it encounters during parsing, print out the internal representation of the modules, and then write it out to a different filename."
        sys.exit(1)
    filename = sys.argv[1]
    
    modules = parse(filename, dump_cleaned = True)

    if True:
        dump_modules(modules)

    convert_to_named_ports(modules)
    convert_to_no_buses(modules)
    convert_to_realistic_fanin(modules)


    write_verilog(modules, filename + ".output.v")



if __name__ == "__main__":
    main()
