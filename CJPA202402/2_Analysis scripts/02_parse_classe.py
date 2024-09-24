import os
import model, ivy
from itertools import chain  # , islice
pd = model.pd
np = pd.np

zero = np.int8(0)
one = np.int8(1)
negone = np.int8(-1)

ALP = model.biomes.states.A

def process_newick(taxon, idx, newick):
    root = ivy.tree.read(newick)
    root.simmap = {}
    root.depth = 0
    for n in root.descendants():
        n.length = np.sum([ x[1] for x in n.simmap ])
        n.depth = n.parent.depth + n.length
    maxdepth = max([ x.depth for x in root.leaves() ])
    for n in root.postiter():
        n.age = maxdepth - n.depth

    d = dict(taxon=taxon, revsample=idx, t=None, node=None, event=0,
             ancestral=0, derived=0, diff=0)

    for a in model.areas:
        d[f'{a.name}'] = zero
        d[f'{a.name}_insitu'] = zero




    for n in root.descendants():
        t = n.parent.age
        # revbayes simmap states are ordered from present to past, so
        # reverse the order to proceed from ancestor to descendant
        simmap = list(reversed(n.simmap))
        s, dur = simmap.pop(0)
        ancidx = int(s)  # index of ancestral state
        ancestral = model.states[ancidx]

        acc = {}
        if n.parent is root:
            # area accumulation
            for a in model.areas:
                acc[f'{a.name}'] = np.int8(
                    bool(ancestral & a) and bool(ancestral & ALP))


        # This is the event for the 'birth' of the branch. Determine
        # if it is an in situ speciation event, i.e.  if `ancestral`
        # is a single area and that area appears in the parent range
        if n.parent and n.parent.simmap:
            # branch begins in a only
            a = (ancestral & model.areas.bitmask)
            # a is in ancestral (parent) range
            b = model.states[int(n.parent.simmap[0][0])]
            if a in model.areas and (a & b):
                name = model.areas(a).name
                if ancestral & ALP:
                    acc[f'{name}'] = one
                    acc[f'{name}_insitu'] = one


        yield {**d, **acc, 't':t, 'node':n.ni, 'event':'birth',
               'ancestral':ancestral, 'derived':ancestral,
               'ancidx':ancidx, 'ancstr':ancestral.name,
               'deridx':ancidx, 'derstr':ancestral.name}
        for s, dur in simmap:
            t -= dur
            deridx = int(s)
            derived = model.states[deridx]
            # change from ancestral to derived at time t
            try:
                e = next(iter(model.stategraph.es.select(
                    _source=ancidx, _target=deridx)))
                event = e['event']
            except StopIteration:
                event = 'cladogenetic'
                daughters = set(
                    chain.from_iterable(model.clado_events_by_idx[ancidx]))
                if deridx not in daughters:
                    event = 'unknown'

            # area accumulation
            acc = {}
            for a in model.areas:
                if (ancestral & a) and not (derived & a):
                    # a is lost
                    if ancestral & ALP:
                        acc[f'{a.name}'] = negone
                elif (not (ancestral & a)) and (derived & a):
                    # a is gained
                    if ancestral & ALP:
                        acc[f'{a.name}'] = one

            # this is the state transition event along the branch
            yield {**d, **acc, 't':t, 'node':n.ni, 'event':event,
                   'ancestral':ancestral, 'derived':derived,
                   'diff':ancestral^derived,
                   'ancidx':ancidx, 'ancstr':ancestral.name,
                   'deridx':deridx, 'derstr':derived.name}
            ancestral = derived
            ancidx = model.states.idx[ancestral]
        # this is the 'death' of the branch
        if not simmap:
            derived = ancestral
            deridx = ancidx
        yield {**d, 't':n.age, 'node':n.ni, 'event':'death',
               'ancestral':ancestral, 'derived':derived,
               'ancidx':ancidx, 'ancstr':ancestral.name,
               'deridx':deridx, 'derstr':derived.name}


def getfile(base, taxon, template):
    for s in taxon, taxon.capitalize(), taxon.lower():
        p = os.path.join(base, taxon, template.format(s))
        if os.path.isfile(p):
            return p
    raise Exception(s)


if __name__ == '__main__':
    ## import ray
    ## ray.init()
    base = 'clades'
    taxa = filter(lambda x: not x.startswith('__'), next(os.walk(base))[1])
    results = {}
    ## @ray.remote
    def process_taxon(taxon):
        Taxon = taxon.capitalize()
        ## tipstates = pd.read_csv(
        ##     getfile(base, taxon, '{}_states.csv'),
        ##     index_col=0, header=None).iloc[:,0]
        simmaps = pd.read_csv(
            getfile(base, taxon, '{}.smap.log'),
            sep='\t', index_col=0).simmap
        N = len(simmaps)
        ## nsamples = 100
        burnin = 0.25
        start = int(N*burnin)
        ## end = N+1
        ## step = int((end-start)/nsamples)
        ## samples = simmaps[start:start+nsamples*step:step]
        samples = simmaps[start:]

        v = []
        for idx, newick in samples.iteritems():
            events = process_newick(Taxon, idx, newick)
            v.append(events)
        df = pd.DataFrame(chain.from_iterable(v))
        print(Taxon, len(v))

        disp_anc = df.ancestral[df.event=='dispersal']
        disp_der = df.derived[df.event=='dispersal']

        disp_anc_nareas = pd.Series(index=disp_anc.index, data=np.sum(
            [ np.bitwise_and(disp_anc, a).astype(bool)
              for a in model.areas ], axis=0))

        # binary columns indicating colonization of H, Y, T ancestors
        for c in 'HBCDEFG':
            area = model.areas.states[c]
            absent_ancestrally = np.logical_not(
                np.bitwise_and(disp_anc, area).astype(bool))
            present_derived = np.bitwise_and(disp_der, area).astype(bool)
            s = np.logical_and(absent_ancestrally, present_derived)
            df[f'{c}_colonization'] = zero
            df[f'{c}_colonization'].update(np.logical_and(
                np.bitwise_and(disp_anc, model.biomes.states.A).astype(bool), s).astype(np.int8))
           
            anc_A = df.ancestral[df[f'{c}_colonization']==1]
            for src in model.areas:
                if src.name == c:
                    continue
                key = f'{c}_colonization_from_{src.name}'
                df[key] = 0.0
                src_in_anc = np.bitwise_and(anc_A, src).astype(bool).astype(np.int8)
                df[key].update(src_in_anc/disp_anc_nareas[anc_A.index])
                

        df.to_csv(f'{Taxon}-events.csv.gz', compression='gzip', index=False)
        del df
        ## return df
        ## results[Taxon] = df

    for taxon in taxa:
        process_taxon(taxon)
