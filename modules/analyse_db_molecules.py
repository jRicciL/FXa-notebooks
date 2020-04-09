from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def load_molecules_from_dir(list_of_sdf_files):
    '''Function to load molecules from sdf files using rdkit'''
    # Load the molecules in a dictionary
    mols_dict = {}
    sanitized = True
    for sdf_file in list_of_sdf_files:
        # Get the molecule name
        mol_name = sdf_file.split('/')[-1].split('.')[0]
        # Try to load the molecule with sanitize = True
        mol_rd = Chem.SDMolSupplier(sdf_file, sanitize = True)[0]
        if mol_rd is None:
            mol_rd = Chem.SDMolSupplier(sdf_file, sanitize = False)[0]
            mol_rd.UpdatePropertyCache(strict = False)
            sanitized = False
        mols_dict[mol_name] = [mol_rd, sanitized]
    return mols_dict

def get_files_list(path_to_sdfs, actives_name = 'ligand', inactives_name = 'decoy'):
    '''Returns a list of path sdf files in a given directory'''
    # Active molecules
    file_list_ACTIVES = glob(path_to_sdfs + F'/{actives_name}*')
    file_list_ACTIVES.sort(key = lambda x: 
                        int(x.split('/')[-1].split('.')[0].split('_')[1]) )
    # Inactive molecules
    file_list_INACTIVES = glob(path_to_sdfs + F'/{inactives_name}*')
    file_list_INACTIVES.sort(key = lambda x: 
                        int(x.split('/')[-1].split('.')[0].split('_')[1]) )
    # Join both list
    file_list = file_list_ACTIVES + file_list_INACTIVES
    return file_list

def get_mol_dataframe(mol_dictionary):
    '''Turns a dictionary of molecules into a dataframe'''
    # Convert to a dataframe
    df = pd.DataFrame(mol_dictionary).T
    df.columns = ['mol_rdk', 'sanitized']
    # Activity inactivity column
    act_inact = ['active' if i[:6] == 'ligand' else 'inactive' for i in df.index]
    df['Activity'] = act_inact
    # Naming the columns
    df = df[['Activity', 'mol_rdk', 'sanitized']]
    return df
