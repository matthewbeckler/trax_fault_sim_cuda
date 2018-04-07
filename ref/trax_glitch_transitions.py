# this is some of the data regarding the potentially glitchy transitions for various gates up to size 8
glitchy_transitions = {}
def dict_multi_insert(dict, keys, value):
    for key in keys:
        #print key, len(value)
        dict[key] = value

#dict_multi_insert(glitchy_transitions, ("and2", "nand2"), [('01', '10'), ('10', '01')]) # 2 transitions
#dict_multi_insert(glitchy_transitions, ("xor2", "xnor2"), [('00', '11'), ('01', '10'), ('10', '01'), ('11', '00')]) # 4 transitions
#dict_multi_insert(glitchy_transitions, ("or2", "nor2"), [('01', '10'), ('10', '01')]) # 2 transitions

# experiments with the new hazard thing to remove the bugs
# basically, when we have transitions like 01->H0 (for and/nand) the 0 masks the H, but not really
dict_multi_insert(glitchy_transitions, ("and2", "nand2"), [('01', '10'), ('10', '01'),
                                                           ('01', 'H0'), ('10', '0H')])
dict_multi_insert(glitchy_transitions, ("or2", "nor2"),   [('01', '10'), ('10', '01'),
                                                           ('01', '1H'), ('10', 'H1')])
# XOR of any X values = X
# XOR of any 0,1,H values = H
# so we don't need to consider H generation here
dict_multi_insert(glitchy_transitions, ("xor2", "xnor2"), [('00', '11'), ('01', '10'), ('10', '01'), ('11', '00')]) # 4 transitions
