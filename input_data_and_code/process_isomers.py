import os
import time
import rdkit
import datetime
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from joblib import Parallel, delayed
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

patterns_to_avoid = ["[R]=[R]=[R]", "[R]#[R]", "[OX2,OX1-]~[OX2,OX1-]"] # array of SMARTS patterns you want to include or exclude in your final data set

def filter_smiles_for_patterns_to_avoid(query_smi, patterns_to_avoid): # function that checks if a given molecule in SMILES format contains any of the patterns previously defined
#                                                                        returns False if so, returns True if no patterns to avoid are present
    try:
        mol = Chem.MolFromSmiles(query_smi)
        for pattern in patterns_to_avoid:
            if mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
                return False
        return True
    except:
        print("Could not run", smi, "through pattern filtering function.")
        
def filter_by_hba_hbd_and_rotors(reference_hba, reference_hbd, reference_rotors, query_smi, margin): # filters molecules by numbers of hydrogen bond acceptors (HBA), donors (HBD), and rotatable bonds with
#                                                                                                      (rotors) with respect to a given reference and margin of tolerance
#                                                                                                      returns None, None, None, None if HBA, HBD or rotors are outside of the given margin, otherwise
#                                                                                                      returns an RDKit molecule object, and the numer of HBA, HBD, and rotors for the query molecule
    try:       
        query_mol = Chem.MolFromSmiles(query_smi)
        query_hba = rdMolDescriptors.CalcNumHBA(query_mol)
        query_hbd = rdMolDescriptors.CalcNumHBD(query_mol)
        query_rotors = rdMolDescriptors.CalcNumRotatableBonds(query_mol)

        if (query_hba <= reference_hba + margin and query_hba >= reference_hba - margin)\
        and (query_hbd <= reference_hbd + margin and query_hbd >= reference_hbd - margin)\
        and (query_rotors <= reference_rotors + margin and query_rotors >= reference_rotors - margin):
            return query_mol, query_hba, query_hbd, query_rotors 
        else:
            return None, None, None, None
    except:
        print("Could not run", query_smi, "through HBA, HBD, and rotors filtering function.")

def organize_data_in_data_frame(i, reference_smi, reference_hba, reference_hbd, reference_rotors, query_smi, query_hba, query_hbd, query_rotors, query_mol, fp1, row, isomers, isomers_from_loose_search):
# function that takes in the numer of the active molecule from actives_final.ism (countung starts from 1), its number of HBA, HBD, and rotors, a query SMILES, its number of HBA, HBD, rotors, an
# RDKit molecule object for the query SMILES string, an RDKit molecular fingeprint object for the active molecule, and two isomers arrays (one containing molecules with an HBA, HBD, and rotors 
# margin of one, the other containing molecules that match a margin of 3. The function computes the Tanimoto similarity between the active molecule and the query and organizes all of the data 
# into a pandas dataframe that gets returned 
    temp = pd.DataFrame(index=range(1), columns=['active', 'active_smi', 'active_HBA', 'active_HBD', 'active_rotors', 'isomer_smi', 'Tanimoto_similarity_to_active', 'isomer_HBA', 'isomer_HBD', 'isomer_rotors'])
    temp.active.iloc[0] = i
    temp.active_smi.iloc[0] = reference_smi
    temp.active_HBA.iloc[0] = reference_hba
    temp.active_HBD.iloc[0] = reference_hbd
    temp.active_rotors.iloc[0] = reference_rotors

    temp.isomer_smi.iloc[0] = query_smi
    temp.isomer_HBA.iloc[0] = query_hba
    temp.isomer_HBD.iloc[0] = query_hbd
    temp.isomer_rotors.iloc[0] = query_rotors

    fp2 = Chem.RDKFingerprint(query_mol)
    tanimoto = DataStructs.TanimotoSimilarity(fp1,fp2)
    temp.Tanimoto_similarity_to_active.iloc[0] = tanimoto

    if row % 100_000 == 0:
        print("Processed", row, "molecules, found", str(len(isomers)) + " suitable isomers from strict search and " + str(len(isomers_from_loose_search)), "suitable isomers from loose search for", i, datetime.datetime.now())
        isomers.to_csv(str(i) + '_isomers.csv', index = False) # save isomers to csv every 100,000 molecules

    return temp
    
# read data for DUD-E actives from actives_final.ism (as distributed by DUD-E) and the corresponding molecular formulae for every molecule in the active set

actives = pd.read_table('actives_final.ism', header = None, sep = '\s+')
formulae = pd.read_table('formulae.txt', header = None, sep = '\s+')

actives.drop(columns = 1, inplace = True)
actives.rename(columns={0:"smiles", 2:"CHEMBL"}, inplace = True)
formulae.rename(columns={1:"formula"}, inplace = True)
formulae.drop(columns=0, inplace = True)
actives = pd.concat([actives, formulae], axis = 1) # group DUD-E active data and molecular formulae into a single pandas data frame
actives.index += 1                                 # set indices to start from 1

def find_suitable_isomers(i, actives, patterns_to_avoid):  # main workhorse function; takes in the number of the active molecule from actives_final.ism (1-based counting), the actives data frame 
#                                                            and any patterns to avoid.  
    if not os.path.isfile('../' + actives.formula.loc[i] + '.smi'): # check if a file of SMILES isomers for the actives exists, return if such a file doesn't exist
        print("No SMILES isomers file for active", i)
        return
        
    print("Now processing active", i, "in", os.getcwd(), datetime.datetime.now())
    row = 0

    if os.path.isfile(str(i) + '.chk'): # check if a checkpoint and csv file from a previous run exist. If so, read in the isomers from the previous run and the number of molecules that have been processed
        f = open(str(i) + '.chk')
        checkpoint = int(f.readline())
        f.close()
        isomers = pd.read_csv(str(i) + '_isomers.csv')
    else:
        checkpoint = 0
        
    isomers = pd.DataFrame(index=range(0), columns=['active', 'active_smi', 'active_HBA', 'active_HBD', 'active_rotors', 'isomer_smi', 'Tanimoto_similarity_to_active', 'isomer_HBA', 'isomer_HBD', 'isomer_rotors'])

    isomers_from_loose_search = pd.DataFrame(index=range(0), columns=['active', 'active_smi', 'active_HBA', 'active_HBD', 'active_rotors', 'isomer_smi', 'Tanimoto_similarity_to_active', 'isomer_HBA', 'isomer_HBD', 'isomer_rotors'])

# calculate HBA, HBD, and rotors for active molecule i (1-based counting) from actives_final.ism and generate an RDKit fingerprint object for it

    reference_smi = actives.smiles.loc[i]
    reference_mol = Chem.MolFromSmiles(reference_smi)
    reference_hba = rdMolDescriptors.CalcNumHBA(reference_mol)
    reference_hbd = rdMolDescriptors.CalcNumHBD(reference_mol)
    reference_rotors = rdMolDescriptors.CalcNumRotatableBonds(reference_mol)       
    fp1 = Chem.RDKFingerprint(reference_mol)    
        
    start = time.time()

    with open('../' + actives.formula.loc[i] + '.smi', 'r') as file: # open MAYGEN-generated isomers file and read in its molecules line-by-line
        for line in file.readlines():
            if (time.time() - start < 7200) and line != '': # a time-limit of 7200 seconds (2 hours) is set for every active. Adjust accordingly. 
                if row < checkpoint:                        # if reading in a molecule that has previously been processed, do nothing 
                    row += 1
                else:
                    try:                                    # if reading in a molecule that has not been processed before, pass it through the function that filters molecules for patterns to avoid
                        query_smi = line.strip()
                        if filter_smiles_for_patterns_to_avoid(query_smi, patterns_to_avoid):
                            query_mol, query_hba, query_hbd, query_rotors = filter_by_hba_hbd_and_rotors(reference_hba, reference_hbd, reference_rotors, query_smi, margin = 1) # if it passes, try filtering for HBA, HBD, and rotors with a margin of 1
                            if (query_mol is not None) and (query_hba is not None) and (query_hbd is not None) and (query_rotors is not None):
                                temp = organize_data_in_data_frame(i, reference_smi, reference_hba, reference_hbd, reference_rotors, query_smi, query_hba, query_hbd, query_rotors, query_mol, fp1, row, isomers, isomers_from_loose_search)
                                isomers = pd.concat([isomers, temp], axis = 0)
                            else:
                                query_mol, query_hba, query_hbd, query_rotors = filter_by_hba_hbd_and_rotors(reference_hba, reference_hbd, reference_rotors, query_smi, margin = 3)  # if it doesn't pass filtering for HBA, HBD, and rotors with a margin of 1, increase the margin to 3 and try again 
                                if (query_mol is not None) and (query_hba is not None) and (query_hbd is not None) and (query_rotors is not None):
                                    temp = organize_data_in_data_frame(i, reference_smi, reference_hba, reference_hbd, reference_rotors, query_smi, query_hba, query_hbd, query_rotors, query_mol, fp1, row, isomers, isomers_from_loose_search)
                                    isomers_from_loose_search = pd.concat([isomers_from_loose_search, temp], axis = 0)

                    except:
                        print("Could not assign RDKit properties for", query_smi, i)
                    row += 1 
            else:
                print("Search ended for", i, "found ", len(isomers), " suitable isomers from strict search and", len(isomers_from_loose_search), "suitable isomers from loose search", datetime.datetime.now())
                break
    
    isomers.drop_duplicates(subset = 'isomer_smi', inplace = True)
    isomers_from_loose_search.drop_duplicates(subset = 'isomer_smi', inplace = True) # explicitly drop duplicates in case there are any
    os.system('echo ' + str(row) + ' > ' + str(i) + '.chk') # save checkpoint file to save the line number, i.e. the number of molecules that have been processed

    if len(isomers) < 100:
        isomers = pd.concat([isomers, isomers_from_loose_search], axis = 0) # if less than 100 molecules with a margin of 1 have been found, add in molecules with a margin of 3.

    isomers.drop_duplicates(subset = 'isomer_smi', inplace = True) # drop duplicates again, just in case there are any
    isomers.sort_values(by = 'Tanimoto_similarity_to_active', ascending = False, inplace = True) # sort by Tanimoto similarity
    isomers.to_csv(str(i) + '_isomers.csv', index = False)                                       # save isomers in csv format

run = Parallel(n_jobs=48)(delayed(find_suitable_isomers)(i, actives, patterns_to_avoid) for i in range(1,len(actives) + 1))  # run using 48 threads, adjust according to available resources

