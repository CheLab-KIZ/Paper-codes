import operator
from functools import reduce
import model

pd = model.pd
areas = model.areas
ranges = model.ranges
biomes = model.biomes
states = model.states
nstates = model.nstates

# read the data csv files - the first column (taxa) is identical in each,
# so we can concatenate them
csvdata = pd.read_csv('range.csv', index_col=0)
csvdata['A']=1

csvdata.iloc[:,:len(areas)] *= areas.values
csvdata.iloc[:,len(areas):] *= biomes.bitvalues
data = pd.Series(
    index=csvdata.index,
    data=[ states.idx[x] for x in csvdata.sum(axis=1).values ])
#Below, change file name and totol species diversity accordingly.
treefile = 'Lycodon.tre' 
outfilebase = 'Lycodon'
total_diversity = 72
nTimeSlices = 1000
ngen_burnin = 200
ngen_mcmc = 100000
revscript = 'Lycodon_classe.rev'

classe_state_datafile = 'Lycodon_states.csv'
data.to_csv(classe_state_datafile, index_label=None)

extinction_rate_lines = []
# extinction rates are nonzero only for single-area states
for v in model.stategraph.vs:
    if v['nareas'] > 1:
        line = f'extinction_rates[{v.index+1}] <- 0.0  # {v["label"]}'
    else:
        line = f'extinction_rates[{v.index+1}] := extinction  # {v["label"]}'
    extinction_rate_lines.append(line)
extinction_rate_lines = '\n'.join(extinction_rate_lines)

# populate non-zero anagenetic rates
anagenetic_rates_lines = []
for e in model.stategraph.es:
    etype = e['event']
    if etype == 'dispersal':
        nsrc = float(len(e['src']))
        line = (f'anarates[{e.source+1}][{e.target+1}] '
                f':= dispersal * {nsrc}  '
                f'# {states[e.source].name}->{states[e.target].name}')
    elif etype == 'extinction':
        line = (f'anarates[{e.source+1}][{e.target+1}] '
                f':= {etype}  '
                f'# {states[e.source].name}->{states[e.target].name}')
    elif etype.startswith('biome'):
        line = (f'anarates[{e.source+1}][{e.target+1}] '
                ':= biome_transition  '
                f'# {states[e.source].name}->{states[e.target].name}')
    anagenetic_rates_lines.append(line)
anagenetic_rates_lines = '\n'.join(anagenetic_rates_lines)

# cladogenetic events
# for each range, enumerate the splits
clado_events_lines = []
extend = clado_events_lines.extend
clado_event_idx = 1
for i, state in enumerate(states):
    # i += 1
    av = list(ranges.decompose(state))  # areas in state
    if len(av) == 1:  # single area range
        extend([f'clado_events[{clado_event_idx}] = [{i},{i},{i}] #[{states[i].name},{states[i].name},{states[i].name}]',
                f'speciation_rates[{clado_event_idx}] := speciation_sym'])
        clado_event_idx += 1
    elif len(av) > 1:
        for a in av:  # iterate over areas in state
            # make single-area state with same biome
            single = a|(state & biomes.bitmask)
            j = states.idx[single]
            clado_events_lines.extend([
                f'clado_events[{clado_event_idx}] = [{i},{i},{j}] #[{states[i].name},{states[i].name},{states[j].name}]',
                f'clado_events[{clado_event_idx+1}] = [{i},{j},{i}] #[{states[i].name},{states[j].name},{states[i].name}]',
                (f'speciation_rates[{clado_event_idx}] '
                 ':= speciation_sub'),
                (f'speciation_rates[{clado_event_idx+1}] '
                 ':= speciation_sub')
                ])
            clado_event_idx += 2
            k = states.idx[state-a]  # index of state without area a
            clado_events_lines.extend([
                f'clado_events[{clado_event_idx}] = [{i},{j},{k}] #[{states[i].name},{states[j].name},{states[k].name}]',
                f'clado_events[{clado_event_idx+1}] = [{i},{k},{j}] #[{states[i].name},{states[k].name},{states[j].name}]',
                (f'speciation_rates[{clado_event_idx}] '
                 ':= speciation_allo'),
                (f'speciation_rates[{clado_event_idx+1}] '
                 ':= speciation_allo')
                ])
            clado_event_idx += 2
clado_events_lines = '\n'.join(clado_events_lines)

unobserved_areas = reduce(operator.ior,
                          [ x for x in areas if csvdata[x.name].sum()==0 ],
                          0)

# root state frequencies
# root_simplex_params = [0]+[1]*(nranges-1)
root_simplex_params = [0]
for i in range(1, len(states)):
    s = states[i]
    if s & unobserved_areas:
        f = 0
    else:
        f = 1.0/len(list(ranges.decompose(s)))
    root_simplex_params.append(f)

with open(revscript, 'w') as outf:
    outf.write(eval(open('classe-template.rev').read()))
