# Helper functions for 2_Docking_analysis folder
import pandas as pd
import numpy as np
from rdkit import Chem
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white', context='talk', font_scale=0.8)



def violin_plot_helper(feature, lig_datasets, xlabel='', ylabel='', title='', figsize=(12,4),
                       split_by_activity=False, palette="Spectral", linewidth=1.8, **kwargs):
    df_ = pd.DataFrame()
    # Create the dataset
    names_ = []
    for name, dataset in lig_datasets.items():
        a = dataset[feature]
        activity = dataset['Activity']
        n_actives = np.sum(activity == 'active')
        length = len(a)
        std_ = np.std(a).round(2)
        mean_ = np.mean(a).round(2)

        names_.append(f'{name}\nn_a/N = {n_actives}/{length}\nMean = {mean_}\nStd = {std_}')

        df_ = df_.append(
                pd.DataFrame(
                    list(zip([name]*length, a, activity)),
                    columns = ['Database', 'Feature', 'Activity']))

    plt.figure(figsize=figsize)
    if split_by_activity:
        _ = sns.violinplot(x='Database', y='Feature', hue = 'Activity', linewidth=linewidth,
                           data=df_, palette=palette, bw=.15, split=split_by_activity, **kwargs)
    else:
        _ = sns.violinplot(x='Database', y='Feature', linewidth=linewidth,
                           data=df_, palette=palette, bw=.15, **kwargs) 

    # plotting
    plt.xticks(np.arange(len(names_)), labels=names_)
    plt.ylabel(ylabel, weight='bold')
    plt.xlabel(xlabel, weight='bold')
    plt.title(title, weight='bold')
    plt.grid(c='lightgrey')

def swarm_plot_helper(feature, lig_datasets, xlabel='', ylabel='', title='', figsize=(12,4),
                       split_by_activity=False, palette="Spectral", linewidth=1.8, **kwargs):
    df_ = pd.DataFrame()
    # Create the dataset
    names_ = []
    for name, dataset in lig_datasets.items():
        a = dataset[feature]
        activity = dataset['Activity']
        n_actives = np.sum(activity == 'active')
        length = len(a)
        std_ = np.std(a).round(2)
        mean_ = np.mean(a).round(2)

        names_.append(f'{name}\nn_a/N = {n_actives}/{length}\nMean = {mean_}\nStd = {std_}')

        df_ = df_.append(
                pd.DataFrame(
                    list(zip([name]*length, a, activity)),
                    columns = ['Database', 'Feature', 'Activity']))

    plt.figure(figsize=figsize)
    if split_by_activity:
        _ = sns.swarmplot(x='Database', y='Feature', hue = 'Activity',
                           data=df_, palette=palette, **kwargs)
    else:
        _ = sns.swarmplot(x='Database', y='Feature', 
                           data=df_, palette=palette,  **kwargs) 

    # plotting
    plt.xticks(np.arange(len(names_)), labels=names_)
    plt.ylabel(ylabel, weight='bold')
    plt.xlabel(xlabel, weight='bold')
    plt.title(title, weight='bold')
    plt.grid(c='lightgrey')

# http://rdkit.blogspot.com/2013/10/comparing-fingerprints-to-each-other.html
'''
def directCompare(scoredLists,fp1,fp2,plotIt=True,silent=False):
    """ Returns: Kendall tau, Spearman rho, and Pearson R
    
    """
    l1 = scoredLists[fp1]
    l2 = scoredLists[fp2]
    rl1=[x[-1] for x in l1]
    rl2=[x[-1] for x in l2]
    vl1=[x[0] for x in l1]
    vl2=[x[0] for x in l2]
    if plotIt:
        _=scatter(vl1,vl2,edgecolors='none')
        maxv=max(max(vl1),max(vl2))
        minv=min(min(vl1),min(vl2))
        _=plot((minv,maxv),(minv,maxv),color='k',linestyle='-')
        xlabel(fp1)
        ylabel(fp2)
    
    tau,tau_p=stats.kendalltau(vl1,vl2)
    spearman_rho,spearman_p=stats.spearmanr(vl1,vl2)
    pearson_r,pearson_p = stats.pearsonr(vl1,vl2)
    if not silent:
        print fp1,fp2,tau,tau_p,spearman_rho,spearman_p,pearson_r,pearson_p
    return tau,spearman_rho,pearson_r'''


from itertools import combinations
from rdkit.DataStructs import FingerprintSimilarity
from rdkit import DataStructs

def compare_lig_db(fp, lig_datasets, method = 'tanimoto', same = None, same_db = ''):
    '''
        Compares pairwise similarity between molecules from two given sets.
    '''
    matched_ligands = {}
    if same:
        combs = (same_db, same_db)
    else:
        combs = combinations(lig_datasets.keys(), 2)
    
    for key_i, key_j in combs:
        print('\n' + '='*20)
        print(key_i, '\t', key_j)
        print('='*20)
        d_i = lig_datasets[key_i]
        d_j = lig_datasets[key_j]

        # Create the list
        matched = []
        for k in d_i.index:
            for p in d_j.index:
                try:
                    fp_sim = FingerprintSimilarity(
                        d_i.loc[k, fp], 
                        d_j.loc[p, fp], metric=DataStructs.TanimotoSimilarity)

                    if fp_sim >= 0.90:
                        # Add to the list
                        matched.append( {'match_mols': (d_i.loc[k, 'mol_rdk'], 
                                                   d_j.loc[p, 'mol_rdk']), 
                                         'match_names': (k, p),
                                         'tanimoto': fp_sim} )
                    if fp_sim >= 0.98:
                        print(k, '\t', p)
                except AttributeError as e:
                    print(e, k, '\t', p)
                    break
        # add to the dict
        matched_ligands[F'{key_i}-{key_j}'] = matched
    return matched_ligands


from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG, Image
from rdkit.Chem import rdDepictor
def draw_matched_ligs(db, matched_ligands, cutoff = 0.99):
    matched_dabatabases = matched_ligands[db]
    mols_to_draw = {}
    for match in matched_dabatabases:
        score = match['tanimoto']
        if score > cutoff:
            # Get both molecules
            mol_i, mol_j = match['match_mols']

            name_i, name_j = match['match_names']
            # Compute 2D coords
            rdDepictor.Compute2DCoords(mol_i)
            rdDepictor.Compute2DCoords(mol_j)

            mols_to_draw[F'{name_i} - {name_j}'] = mol_i
            print('='*25 + ' '*5 + db +  ' '*5 + '='*25)
            img = Chem.Draw.MolsToGridImage(
                (mol_i, mol_j), legends = (name_i, name_j),
                molsPerRow = 2, subImgSize = (300,200))
            display(img)
            mol_i.GetSubstructMatch(mol_j)
            display(mol_i)