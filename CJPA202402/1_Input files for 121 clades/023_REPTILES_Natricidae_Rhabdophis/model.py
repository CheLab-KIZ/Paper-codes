import igraph
import pandas as pd
from bitstates import Bitstates, IntFlag
from collections import defaultdict
np = pd.np

areas = Bitstates('area', list('HBCDEFG'))
ranges = Bitstates('range', areas.bitnames, data='compound', null=True)
biomes = Bitstates('biome', list('A'), data='compound',
                   offset=len(ranges.bitnames))
null = ranges[0]
v = [(null.name, null.value)]
nbits = {null.name:0, 0:0}
for b in biomes:
    for r in ranges[1:]:
        name = f'{r.name}+{b.name}'
        value = b|r
        v.append((name, value))
        nb = ranges.nbits[r]+biomes.nbits[b]
        nbits[name] = nb
        nbits[value] = nb

# here, states are a combination of geographic range and biome
states = Bitstates('state', [ x[0] for x in v ],
                   data=(IntFlag('state', v), nbits))
nstates = len(states)

# in this simple case, all areas are connected
areagraph = igraph.Graph.Full(len(areas), directed=False)
areagraph.vs['bits'] = areas.bitvalues

# make graph connecting states through anagenetic range and
# biome transitions
def connected(g, bits1, bits2):
    '''
    determine if the vertices defined by bits1 and bits2 are connected in
    graph g
    '''
    try:
        g.get_eid(areagraph.vs.select(bits=bits1)[0],
                  areagraph.vs.select(bits=bits2)[0])
        return True
    except:
        return False


g = igraph.Graph(directed=True)
g.add_vertices(len(states))
g.vs['label'] = states.names
g.vs['nareas'] = [ len(tuple(ranges.decompose(x))) for x in states ]
# add edges from single-area states to null (globally extinct)
for a in areas:
    for b in biomes:
        i = states.idx[a|b]
        g.add_edge(i, 0, event='extinction')
for i, src in enumerate(states[:-1]):
    for j, tgt in enumerate(states[i+1:]):
        source = src; target = tgt
        j += i+1
        bit = src^tgt  # the bits that differ
        if bit in areas.values:  # single bit
            a = areas.states(bit)
            si, ti = i, j
            if tgt & bit:  # tgt has extra bit
                source, target = target, source
                si, ti = j, i
            # store bits in src connected to bit
            target = states[ti]
            v = [ x for x in ranges.decompose(target)
                  if connected(areagraph, bit, x) ]
            if v:
                g.add_edge(si, ti, event='extinction', area=areas.states(bit))
                g.add_edge(ti, si, event='dispersal', src=v)
        elif bit in biomes.bitvalues:
            biome = biomes.states(bit)
            si, ti = i, j
            if tgt & bit:  # tgt has extra bit
                source, target = target, source
                si, ti = j, i
            g.add_edge(si, ti, event='biome_contraction')
            g.add_edge(ti, si, event='biome_expansion')
stategraph = g

clado_events_by_idx = defaultdict(set)
clado_events_by_state = defaultdict(set)
for i, state in enumerate(states):
    # i += 1
    av = list(ranges.decompose(state))  # areas in state
    if len(av) == 1:  # single area range
        clado_events_by_idx[i].add(frozenset([i]))
        clado_events_by_state[state].add(frozenset([state]))
    elif len(av) > 1:
        for a in av:  # iterate over areas in state
            # make single-area state with same biome
            single = a|(state & biomes.bitmask)
            j = states.idx[single]
            clado_events_by_idx[i].add(frozenset((i,j)))
            clado_events_by_state[state].add(frozenset((state,single)))
            k = states.idx[state-a]  # index of state without area a
            clado_events_by_idx[i].add(frozenset((j,k)))
            clado_events_by_state[state].add(frozenset((single,state-a)))
