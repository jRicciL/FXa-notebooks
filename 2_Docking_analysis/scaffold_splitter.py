import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold 
from rdkit.Chem import MolToSmiles
from rdkit.Chem import rdDepictor

to2d = Chem.rdDepictor.Compute2DCoords

# Define a lambda function to compute Murcko scaffolds
def scaffold2smiles(mol, generic=True, return_smiles = True):
    ''' Returns a SMILES string representing the Murcko Scaffold of a given molecule'''
    if generic:
        # Makes a Murcko scaffold generic (all atoms -> carbon and all bonds -> single)
        scff = MurckoScaffold.MakeScaffoldGeneric(mol)
        scff = MurckoScaffold.GetScaffoldForMol(scff)
        scff_smiles = MolToSmiles(scff)
    else:
        # Return a smiles scaffold 
        try:
            scff = MurckoScaffold.GetScaffoldForMol(mol)
            scff_smiles =  MolToSmiles(scff)
        except:
            scff_smiles = '' 
            scff = np.nan
    if return_smiles:
        return scff_smiles
    else:
        return scff


def scaffold_splitter(scaffold_series, test_size=0.2, stratify=None):
    '''
    Performs a Train Test splitting using Murcko Scaffolds.
    
    Parameters:
    -----------
        scaffold_series: pandas Series.
            Series object with precomputed Murcko Scaffolds (SMILES format).
        train_size, test_size, valid_size: float.
            Train, Validation and Test fraction size for each subset.
    Retunrs:
    --------
        train_inds, valid_inds, test_inds: array-like.
            Arrays containing the molecule indices of the corresponding train, test subset.
    '''
    
    assert test_size < 0.8, 'test_size must be a float number less than 0.8'
    
    train_size = 1 - test_size
    scaffolds = {}
    train_inds, valid_inds, test_inds = [], [], []
    
    if isinstance(stratify, pd.Series): 
        # Use stratify (y) to get active indices
        actives = (stratify == 1)
        scaffold_series_actives = scaffold_series[actives] # Keep Actives
        scaffold_series = scaffold_series[~actives] # Keep Inactives
        # Perform Train/Test Splitting only for actives
        act_train_inds, act_valid_inds, act_test_inds = scaffold_splitter(
                                              scaffold_series_actives, 
                                              test_size=test_size, stratify=None)
        # Start filling train_inds, valid_inds, test_inds with active indices
        train_inds, valid_inds, test_inds = act_train_inds, act_valid_inds, act_test_inds
    
    data_len = len(scaffold_series)
    
    for ind, scff in scaffold_series.items():
        if scff not in scaffolds:
                scaffolds[scff] = [ind]
        else:
            scaffolds[scff].append(ind)
    # Sort from largest t smallest scaffold sets
    scaffolds = {key: sorted(value) for key, value in scaffolds.items()}
    # Scaffold sets
    scaffold_sets = [
        scff_set
        for (scff, scff_set) in sorted(
            scaffolds.items(), key=lambda x: (len(x[1]), x[1][0]), reverse=True)
    ]
    
    # create Train and Test sets
    train_cutoff = int(train_size * data_len)
    
    # Full Train-Test sets
    for scff_set in scaffold_sets:
        if len(train_inds) + len(scff_set) > train_cutoff:
            test_inds += scff_set
        else:
            train_inds += scff_set
    return train_inds, valid_inds, test_inds
    
    
def train_test_scaffold_split(X, y, scaffold_series, test_size=0.2, stratify=None):
    '''
    Performs a Train Test splitting using Murcko Scaffolds.
    
    Parameters:
    -----------
        X: pandas DataFrame
            Dataframe with m rows and n columns.
        y: pandas Series
            Series object with m target values. Index values should match those in X.
        scaffold_series: pandas Series.
            Series object with precomputed Murcko Scaffolds (SMILES format).
        test_size: float
            Test seti fraction size.
        stratify: None or pandas Series
            Series object with target values to stratify on into Train and Test sets.
        
    Retunrs:
    --------
        X_train, X_test: pandas DataFrame.
            Train and Test X subset accordingly to the Scaffold Splitting.
        y_train, y_test: pandas Series.
            Train and Test target values.
    '''

    assert (X.index == scaffold_series.index).all(), 'Index should be the same between X and rdk_mols_series'
    assert len(X) == len(y), 'Numer of observations in X and y is not the same.'
    
    train_inds, _, test_inds = scaffold_splitter(scaffold_series,
                                    test_size=test_size, stratify=stratify)
    # Subset X and y
    X_train = X.loc[train_inds]
    y_train = y.loc[train_inds]
    X_test = X.loc[test_inds]
    y_test = y.loc[test_inds]
    
    return X_train, X_test, y_train, y_test 