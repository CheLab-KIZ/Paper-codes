from enum import IntFlag
from itertools import combinations
from collections.abc import Sequence
from functools import reduce

class Bitstates(Sequence):
    """
    A class to represent an enumeration of categorical states, each
    stored as a bit using Python's enum.IntFlag.
    """
    def __init__(self, name, bitnames, data='simple', null=False, sep='',
                 maxbits=None, offset=0):
        self.name = name
        self.bitnames = bitnames
        self.null = null
        self.sep = sep
        self.offset = offset
        self.maxbits = maxbits
        if data == 'simple':
            states, nbits = Bitstates.generate_simple(
                name, bitnames, null, offset)
        elif data == 'compound':
            states, nbits = Bitstates.generate_compound(
                name, bitnames, null, sep, maxbits, offset)
        else:  # data = (states, nbits)
            states, nbits = data
        self.states = states
        self.nbits = nbits
        self.list = list(states)
        self.size = len(self.states)
        self.names = tuple([ x.name for x in self.states ])
        self.values = tuple([ x.value for x in self.states ])
        self.idx = {}
        for i, x in enumerate(self.states):
            self.idx[x.name] = i
            self.idx[x.value] = i
        self.bitvalues = tuple([ states[x] for x in bitnames ])
        self.bitmask = 0
        for b in self.bitvalues:
            self.bitmask |= b

    def decompose(self, value):
        for v in self.bitvalues:
            if value & v:
                yield v

    ## def shift(self, offset):
    ##     self.offset += offset
    ##     if offset > 0:
    ##         self.values = [ x << offset for x in self.values ]
    ##     else:
    ##         self.values = [ x >> -offset for x in self.values ]
    ##     self.states = IntFlag(self.name, zip(self.names, self.values))
    ##     self.list = list(self.states)
    ##     for k in list(self.nbits):
    ##         if isinstance(k, int):
    ##             nb = self.nbits[k]
    ##             i = self.idx[k]
    ##             del self.nbits[k]
    ##             del self.idx[k]
    ##             if offset > 0:
    ##                 self.nbits[k << offset] = nb
    ##                 self.idx[k << offset] = i
    ##             else:
    ##                 self.nbits[k >> -offset] = nb
    ##                 self.idx[k >> -offset] = i

    @staticmethod
    def generate_simple(name, bitnames, null=False, offset=0):
        flags = []
        if null:
            flags.append(('null', 0))
        for i, s in enumerate(bitnames):
            i += offset
            flags.append((s, 1 << i))
        states = IntFlag(name, flags)
        nbits = {}
        for state in states:
            nbits[state.name] = 1
            nbits[state.value] = 1
        return states, nbits

    @staticmethod
    def generate_compound(name, bitnames, null=False, sep='', maxbits=None,
                          offset=0):
        for s in bitnames:
            if len(s) > 1:
                sep = '+'
                break
        nbits = {}
        s2i = dict([ (s, i+offset) for i, s in enumerate(bitnames) ])
        d = []
        if null:
            d.append(('null', 0))
            nbits['null'] = 0
            nbits[0] = 0
        for s in bitnames:
            v = 1 << s2i[s]
            d.append((s, v))
            nbits[s] = 1
            nbits[v] = 1
        for i in range(2, len(bitnames)+1):
            if maxbits and i > maxbits:
                continue
            for c in combinations(bitnames, i):
                s = sep.join(c)
                j = 0
                for x in c:
                    j |= 1 << s2i[x]
                d.append((s, j))
                nbits[s] = i
                nbits[j] = i
        states = IntFlag(name, d)
        return states, nbits

    @staticmethod
    def generate_from_graph(name, bitnames, g, null=True, sep='',
                            maxbits=None):
        '''
        calculate the set of "contiguous" compound states with maximum
        size `maxbits` from the given graph of simple states
        '''
        # http://stackoverflow.com/questions
        # /15658245/efficiently-find-all-connected-subgraphs
        for s in bitnames:
            if len(s) > 1:
                sep = '+'
                break
        nbits = {}
        results = []
        if null:
            results.append(('null', 0))
            nbits['null'] = 0
            nbits[0] = 0
        for i,s in enumerate(bitnames):
            j = 1 << i
            results.append((s, j))
            nbits[s] = 1
            nbits[j] = 1

        def gsg(nodes, subset, neighbors):
            if not subset:
                candidates = nodes
            else:
                candidates = nodes.intersection(neighbors)
            if not candidates:
                nb = len(subset)
                if 2 <= nb <= maxbits:
                    s = sep.join([ bitnames[x.index] for x in
                                   sorted(subset, key=lambda e:e.index) ])
                    j = 0
                    for x in subset:
                        j |= 1 << x.index
                    results.append((s, j))
                    nbits[s] = nb
                    nbits[j] = nb
            else:
                n = next(iter(candidates))
                sn = set([n])
                nmsn = nodes - sn
                gsg(nmsn, subset, neighbors)
                gsg(nmsn, subset.union(sn), neighbors.union(n.neighbors()))

        gsg(set(g.vs()), set(), set())
        states = IntFlag(name, results)
        return states, nbits

    def __getitem__(self, i):
        return self.list[i]

    def __iter__(self):
        return iter(self.states)

    def __len__(self):
        return self.size

    def __call__(self, x):
        return self.states(x)
    
