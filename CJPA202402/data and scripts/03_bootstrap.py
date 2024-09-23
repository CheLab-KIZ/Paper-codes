# make a bootstrapped dataframe from the bootstrap index file
import pandas as pd
import os
import glob
np = pd.np
cwd = os.getcwd()

#os.chdir(read_path)
path =r'' # use your path
allFiles = glob.glob(path + "/*events.csv.gz")

frame = pd.DataFrame()
for file_ in allFiles:
    df = pd.read_csv(file_,index_col=None, header=0)
    frame=frame.append(df)
#frame.to_csv('events.csv.gz', compression='gzip', index=False)
frame.set_index(['taxon', 'revsample'], inplace = True,drop = False)
frame.taxon.unique()
frame.sort_index(inplace=True)
uidx = frame.index.drop_duplicates()
taxa = uidx.get_level_values(0).drop_duplicates()

# make a dataframe where taxa are columns and rows are random revsample values,
# i.e. each row is a bootstrap replicate
bootdf = pd.DataFrame(dict(
    (t, np.random.choice(uidx[uidx.get_loc(t)].get_level_values(1), 1000))
    for t in taxa))
bootdf.to_csv('bootstrap-idx.csv', index=False)

#dtype = eval(open('events.dtypes').read())
#df = pd.read_csv('events.csv.gz', index_col=['taxon', 'revsample'], dtype=dtype)
#df.sort_index(inplace=True)
frame.sort_index(inplace=True)
print('sorted events')
bootidx = pd.read_csv('bootstrap-idx.csv')

def pull_rep(bootrep, row):
    # pull out one random stochastic map for each taxon
    sample = frame.loc[row.items()].sort_values('t')
    sample['bootrep'] = bootrep
    return sample

def pull():
    for bootrep, row in bootidx.iterrows():
        print('pulled', bootrep)
        yield pull_rep(bootrep, row)

bootdf = pd.concat(pull(), ignore_index=True)

bootdf.to_csv('events-bootstrap1000.csv.gz',index=False, compression='gzip')
