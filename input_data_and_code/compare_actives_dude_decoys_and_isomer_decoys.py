import os
import re
import time
import rdkit
import datetime
import numpy as np
import pandas as pd
from glob import glob
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from scipy.stats import linregress
from rdkit.Chem import Descriptors
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

folders = glob("*/")

def compare_ligs(folder): # for every target (ada, comt, cxcr4, def, fabp4, fak1, glcm, grik1, hivint, hs90a, hxk4, mapk2, mcr, nos1, nram, pa2ga, pur2, pygm, rock1, sahh) which is a separate folder,
# compute DUD-E actives - DUD-E actives, DUD-E decoys - DUD-E decoys, and isomer decoys - isomer decoys Tanimoto similarities and plot them 

    os.chdir(folder)
    print(os.getcwd())
    
    actives = pd.read_table('actives_final.ism', header = None, sep = '\s+')    # read in data for DUD-E actives and DUD-E decoys
    dude_decoys = pd.read_table('decoys_final.ism', header = None, sep = '\s+')

    actives.rename(columns = {0:'smiles', 2:"CHEMBLID"}, inplace = True)
    actives.drop(columns = 1, inplace = True)

    dude_decoys.rename(columns = {0:'smiles', 1:"ID"}, inplace = True)

    actives.index += 1     # set indices to start from 1
    dude_decoys.index += 1 # set indices to start from 1

    for i in range(1,len(actives) + 1):
        if os.path.isfile(str(i) + '_isomers.csv'):
            data = pd.read_csv(str(i) + '_isomers.csv')
            data.drop_duplicates(subset = 'isomer', inplace = True)
            data.sort_values(by = 'Tanimoto_similarity', ascending = False, inplace = True)
            data.reset_index(inplace = True, drop = True)

        if len(data) < 100:
            actives.drop(index = i, inplace = True) # drop actives that have less than 100 isomers
        else:
            data = data[0:100]

            data.to_csv(str(i) + '_100_isomers.csv', index = False)
            
            dude_actives_dude_actives_tanimotos = []

            dude_decoys_dude_decoys_tanimotos = []
            
            isomer_decoys_isomer_decoys_tanimotos = []

    for active in actives.index: # compute DUD-E active - DUD-E active Tanomotos
            mol = Chem.MolFromSmiles(actives.smiles.loc[active])
            fp_active = Chem.RDKFingerprint(mol)
            
            for another_active in actives.index:
                if active != another_active:
                    another_mol = Chem.MolFromSmiles(actives.smiles.loc[another_active])
                    fp_another_active = Chem.RDKFingerprint(another_mol)
                    dude_actives_dude_actives_tanimotos.append(DataStructs.TanimotoSimilarity(fp_active,fp_another_active)) 
                    
    print("Done computing active - active Tanomotos")                    
    np.save('dude_actives_dude_actives_tanimotos.npy', dude_actives_dude_actives_tanimotos) # save DUD-E active - DUD-E active Tanomotos to .npy file for safe keeping

    for decoy in dude_decoys.index: # compute DUD-E decoy - DUD-E decoy Tanomotos
            mol = Chem.MolFromSmiles(dude_decoys.smiles.loc[decoy])
            fp_decoy = Chem.RDKFingerprint(mol)
            
            for another_decoy in dude_decoys.index:
                if decoy != another_decoy:
                    another_mol = Chem.MolFromSmiles(dude_decoys.smiles.loc[another_decoy])
                    fp_another_decoy = Chem.RDKFingerprint(another_mol)
                    dude_decoys_dude_decoys_tanimotos.append(DataStructs.TanimotoSimilarity(fp_decoy,fp_another_decoy)) 
                    
            print("Done with DUD-E", decoy, "out of", len(dude_decoys), "DUD-E decoys")

    print("Done computing DUD-E decoy - DUD-E decoy Tanomotos")                    
    np.save('dude_decoys_dude_decoys_tanimotos.npy', dude_decoys_dude_decoys_tanimotos) # save DUD-E decoy - DUD-E decoy Tanomotos to .npy file for safe keeping                
                
    csvs = glob("*100_isomers.csv") 
    all_isomer_data = pd.DataFrame(index = range(0), columns = ['isomer','Tanimoto_similarity'])    
    
    for csv in csvs:
        temp = pd.read_csv(csv)
        all_isomer_data = pd.concat([all_isomer_data, temp], axis = 0)

    all_isomer_data.reset_index(inplace = True, drop = True)
    
    for decoy in all_isomer_data.index: # compute isomer decoy - isomer decoy Tanomotos
            mol = Chem.MolFromSmiles(all_isomer_data.isomer.iloc[decoy])
            fp_decoy = Chem.RDKFingerprint(mol)
            
            for another_decoy in all_isomer_data.index:
                if decoy != another_decoy:
                    another_mol = Chem.MolFromSmiles(all_isomer_data.isomer.iloc[another_decoy])
                    fp_another_decoy = Chem.RDKFingerprint(another_mol)
                    isomer_decoys_isomer_decoys_tanimotos.append(DataStructs.TanimotoSimilarity(fp_decoy,fp_another_decoy)) 

            print("Done with isomer", decoy, "out of", len(all_isomer_data), "isomer decoys")

    print("Done computing isomer decoy - isomer decoy Tanomotos")                    
    np.save('isomer_decoys_isomer_decoys_tanimotos.npy', isomer_decoys_isomer_decoys_tanimotos) # save isomer decoy - isomer decoy Tanomotos to .npy file for safe keeping               
                

    
###################################################### PLOT TANIMOTOS FOR INDIVIDUAL TARGETS ########################################################

    plt.title('Tanimoto similarities for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(dude_actives_dude_actives_tanimotos, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - experimental actives similarities')
    plt.hist(dude_decoys_dude_decoys_tanimotos, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental decoys - experimental decoys similarities')
    plt.hist(isomer_decoys_isomer_decoys_tanimotos, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization decoys - isomerization decoys similarities')
    plt.legend(fontsize=11)
    plt.xlabel('Tanimoto similarity', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_distributions.jpeg', dpi = 600)
    plt.close()

    os.chdir('..')

    return dude_actives_dude_actives_tanimotos, dude_decoys_dude_decoys_tanimotos, isomer_decoys_isomer_decoys_tanimotos 
         
run = Parallel(n_jobs=20)(delayed(compare_ligs)(folders[i]) for i in range(0,len(folders))) # run analysis on all 20 targets simultaneously and aggregate the data

dude_actives_dude_actives_tanimotos = list(run[0][0]) + list(run[1][0]) + list(run[2][0]) + list(run[3][0]) + list(run[4][0]) + list(run[5][0]) + list(run[6][0]) + list(run[7][0]) + list(run[8][0]) + list(run[9][0]) + list(run[10][0]) + list(run[11][0]) + list(run[12][0]) + list(run[13][0]) + list(run[14][0]) + list(run[15][0]) + list(run[16][0]) + list(run[17][0]) + list(run[18][0]) + list(run[19][0])

dude_decoys_dude_decoys_tanimotos = list(run[0][1]) + list(run[1][1]) + list(run[2][1]) + list(run[3][1]) + list(run[4][1]) + list(run[5][1]) + list(run[6][1]) + list(run[7][1]) + list(run[8][1]) + list(run[9][1]) + list(run[10][1]) + list(run[11][1]) + list(run[12][1]) + list(run[13][1]) + list(run[14][1]) + list(run[15][1]) + list(run[16][1]) + list(run[17][1]) + list(run[18][1]) + list(run[19][1])

isomer_decoys_isomer_decoys_tanimotos = list(run[0][2]) + list(run[1][2]) + list(run[2][2]) + list(run[3][2]) + list(run[4][2]) + list(run[5][2]) + list(run[6][2]) + list(run[7][2]) + list(run[8][2]) + list(run[9][2]) + list(run[10][2]) + list(run[11][2]) + list(run[12][2]) + list(run[13][2]) + list(run[14][2]) + list(run[15][2]) + list(run[16][2]) + list(run[17][2]) + list(run[18][2]) + list(run[19][2])

###################################################### PLOT TANIMOTOS FOR ENTIRE DATA SET ########################################################

plt.title('Tanimoto similarities', fontsize = 20)
plt.hist(dude_actives_dude_actives_tanimotos, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - experimental actives similarities')
plt.hist(dude_decoys_dude_decoys_tanimotos, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental decoys - experimental decoys similarities')
plt.hist(isomer_decoys_isomer_decoys_tanimotos, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization decoys - isomerization decoys similarities')
plt.legend(fontsize=11)
plt.xlabel('Tanimoto similarity', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_distributions.jpeg', dpi = 600)
plt.close()


