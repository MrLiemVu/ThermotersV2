import numpy as np
from collections import OrderedDict

bases = "acgt"
lett2idx = dict(zip(bases, range(4)))

def getRBSpositions(data):
    '''
    Find the position of the RBS in the sequence.
    
    Inputs:
        data: pandas DataFrame
        
    Outputs:
        rbspos: int
    '''
    rbs_seq = "aggag" #"aggagaag"
    seq = data.iloc[0, 0]
    rbspos = seq.find(rbs_seq)
    return rbspos

def splitDataset(data,
    split_fractions = np.array([.6, .2, .2]),
    split_names = ["training", "validation", "evaluation"],
    rndSeed = 1):
    '''
    Split data into training, validation and evaluation sets.
    
    Inputs:
        data: pandas DataFrame
        split_fractions: array of floats summing to 1
        split_names: list of strings
        rndSeed: int
    
    Outputs:
        data: pandas DataFrame with new column "data split"
    '''
    assert split_fractions.sum() == 1
    data["data split"] = ""
    Ndata = len(data)
    cut_indices = np.cumsum(split_fractions * Ndata)[:-1].astype(int)
    np.random.seed(rndSeed)
    rnd_idx = np.random.permutation(data.index)
    splits = np.split(rnd_idx, cut_indices)
    for name, idx in zip(split_names, splits):
        data.loc[idx,"data split"] = name

    return data

def numerizeSequences(sequences, Ltot=105, pad="a", start=0):
    '''
    Convert a list of sequences to a numerical representation.
    
    Inputs:
        sequences: list of strings
        Ltot: int  # Total length of the sequence
        pad: string  # Padding character
        start: int  # Start position
        
    Outputs:
        out: numpy array of integers
    '''
    from pandas import Series

    out = np.zeros((len(sequences), Ltot), dtype=np.int8)
    lengths = list(map(len, sequences))
    mostCommonLength = Series(lengths).value_counts().idxmax()

    for i, seq in enumerate(sequences):
        if len(seq)<mostCommonLength:
            seq = (mostCommonLength-len(seq)) * pad + seq
        out[i] = [lett2idx[l] for l in seq[start:start+Ltot]]

    return out

def createNumData(
    DataDict,
    tts=["training", "validation", "evaluation"],
    rndSeed=1
):
    '''
    Create numerical data from a dictionary of dataframes.

    Inputs:
        DataDict: dictionary of pandas DataFrames
        tts: list of strings
        rndSeed: int
    
    Outputs:
        numDataDict: dictionary of dictionaries of numpy arrays
    '''
    numDataDict = OrderedDict([(k, OrderedDict()) for k in tts])
    for dataID in DataDict:
        expressions = DataDict[dataID].iloc[:,1]
        Nbins = int(np.round(expressions).max() +1 )
        for tt in tts:
            if tt == "all":
                data = DataDict[dataID]
            else:
                fltr = DataDict[dataID]["data split"].isin(tt.split("+"))
                data = DataDict[dataID][fltr]
            seqs = numerizeSequences(data['sequence'],
#                                      rcrop = 10 if dataID=="36N" else 14,
#                                      Ltot = 100 if dataID=="36N" else 120
                                     Ltot = 115
                                    )
            lums = data.iloc[:, 1]
            np.random.seed(rndSeed)
            jitter = (np.random.rand(len(lums)) - .5) * 1e-9
            digiLums = np.round(lums + jitter).astype(int)
            weights = np.nan * np.ones_like(lums)
            for jb in range(Nbins):
                fltr = digiLums == jb
                if sum(fltr):
                    weights[fltr] = 1.0 / sum(fltr)
            weights /= weights.sum()
            weights *= len(weights)
            numDataDict[tt][dataID] = {
                "seqs":seqs,
                "lums":lums,
                "digiLums":digiLums,
                "weights":weights
            }
    return numDataDict