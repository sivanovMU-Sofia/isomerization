import os
import time
import rdkit
import datetime
import numpy as np
import pandas as pd
from glob import glob
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

from scipy.stats import mannwhitneyu

folders = glob("*/")

def compare_ligs(folder): # for every target (ada, comt, cxcr4, def, fabp4, fak1, glcm, grik1, hivint, hs90a, hxk4, mapk2, mcr, nos1, nram, pa2ga, pur2, pygm, rock1, sahh) which is a separate folder,
# compute molecular properties for DUD-E actives, DUD-E decoys, and isomer decoys, and plot them 
    os.chdir(folder)
    print(os.getcwd())
    
    actives = pd.read_table('actives_final.ism', header = None, sep = '\s+')    # read in DUD-E actives
    dude_decoys = pd.read_table('decoys_final.ism', header = None, sep = '\s+') # read in DUD-E decoys

    actives.rename(columns = {0:'smiles', 2:"CHEMBLID"}, inplace = True)
    actives.drop(columns = 1, inplace = True)

    dude_decoys.rename(columns = {0:'smiles', 1:"ID"}, inplace = True)

    actives.index += 1      # set indices to start from 1
    dude_decoys.index += 1  # set indices to start from 1

    actives['MolWt'] = np.nan
    actives['HBA'] = np.nan
    actives['HBD'] = np.nan
    actives['rotors'] = np.nan
    actives['charge'] = np.nan
    actives['TPSA'] = np.nan
    actives['LabuteASA'] = np.nan

    dude_decoys['MolWt'] = np.nan
    dude_decoys['HBA'] = np.nan
    dude_decoys['HBD'] = np.nan
    dude_decoys['rotors'] = np.nan
    dude_decoys['charge'] = np.nan
    dude_decoys['TPSA'] = np.nan
    dude_decoys['LabuteASA'] = np.nan

    for j in actives.index: # compute molecular weights, HBA, HBD, rotors, charges, TPSA, and LabuteASA for DUD-E actives
        mol = Chem.MolFromSmiles(actives.smiles.loc[j])
        actives['MolWt'].loc[j] = Descriptors.MolWt(mol)
        actives['HBA'].loc[j] = rdMolDescriptors.CalcNumHBA(mol)
        actives['HBD'].loc[j] = rdMolDescriptors.CalcNumHBD(mol)
        actives['rotors'].loc[j] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        AllChem.ComputeGasteigerCharges(mol)
        charge = 0
        for atom in mol.GetAtoms():
            charge += atom.GetDoubleProp("_GasteigerCharge")
        actives['charge'].loc[j] = charge
        actives['TPSA'].loc[j] = rdMolDescriptors.CalcTPSA(mol)
        actives['LabuteASA'].loc[j] = rdMolDescriptors.CalcLabuteASA(mol)

    for j in dude_decoys.index: # compute molecular weights, HBA, HBD, rotors, charges, TPSA, and LabuteASA for DUD-E decoys
        mol = Chem.MolFromSmiles(dude_decoys.smiles.loc[j])
        dude_decoys['MolWt'].loc[j] = Descriptors.MolWt(mol)
        dude_decoys['HBA'].loc[j] = rdMolDescriptors.CalcNumHBA(mol)
        dude_decoys['HBD'].loc[j] = rdMolDescriptors.CalcNumHBD(mol)
        dude_decoys['rotors'].loc[j] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        AllChem.ComputeGasteigerCharges(mol)
        charge = 0
        for atom in mol.GetAtoms():
            charge += atom.GetDoubleProp("_GasteigerCharge")
        dude_decoys['charge'].loc[j] = charge
        dude_decoys['TPSA'].loc[j] = rdMolDescriptors.CalcTPSA(mol)
        dude_decoys['LabuteASA'].loc[j] = rdMolDescriptors.CalcLabuteASA(mol)

    for i in range(1,len(actives) + 1): # read in data for MAYGEN isomers
        if os.path.isfile(str(i) + '_isomers.csv'):
            data = pd.read_csv(str(i) + '_isomers.csv')
            data.drop_duplicates(subset = 'isomer_smi', inplace = True)
            data.reset_index(inplace = True, drop = True)

        if len(data) < 100: # if less than 100 isomers for a given active are present, remove the active from consideration, else take the top 100 (by Tanimoto similarity) actives
            actives.drop(index = i, inplace = True)
        else:
            data = data[0:100]
                
            mol = Chem.MolFromSmiles(data.active_smi.iloc[0]) # for the 100 actives, compute molecular weight, HBA, HBD, rotors, charges, TPSA, and LabuteASA
            data['active_MolWt'] = Descriptors.MolWt(mol)
            AllChem.ComputeGasteigerCharges(mol)
            charge = 0
            for atom in mol.GetAtoms():
                charge += atom.GetDoubleProp("_GasteigerCharge")
            data['active_charge'] = charge
            data['active_TPSA'] = rdMolDescriptors.CalcTPSA(mol)
            data['active_LabuteASA'] = rdMolDescriptors.CalcLabuteASA(mol)

            data['isomer_MolWt'] = np.nan
            data['isomer_charge'] = np.nan
            data['isomer_TPSA'] = np.nan
            data['isomer_LabuteASA'] = np.nan

            data = data[['active', 'active_smi', 'active_MolWt', 'active_HBA', 'active_HBD', 'active_rotors', 'active_charge', 'active_TPSA', 'active_LabuteASA', 'isomer_smi', 'Tanimoto_similarity_to_active', 'isomer_MolWt', 'isomer_HBA', 'isomer_HBD', 'isomer_rotors', 'isomer_charge', 'isomer_TPSA', 'isomer_LabuteASA']]

            for j in data.index:
                mol = Chem.MolFromSmiles(data.isomer_smi.iloc[j])
                data.isomer_MolWt.iloc[j] = Descriptors.MolWt(mol)
                AllChem.ComputeGasteigerCharges(mol)
                charge = 0
                for atom in mol.GetAtoms():
                    charge += atom.GetDoubleProp("_GasteigerCharge")
                data.isomer_charge.iloc[j] = charge
                data.isomer_TPSA.iloc[j] = rdMolDescriptors.CalcTPSA(mol)
                data.isomer_LabuteASA.iloc[j] = rdMolDescriptors.CalcLabuteASA(mol)

            data.to_csv(str(i) + '_100_isomers.csv', index = False) # save data for the 100 isomers for every active to a csv file

            dude_actives_dude_decoys_tanimotos = []


    for active in actives.index: # compute Tanimoto similarities for every DUD-E active - DUD-E decoy pair
            mol = Chem.MolFromSmiles(actives.smiles.loc[active])
            fp_active = Chem.RDKFingerprint(mol)
            
            for decoy in dude_decoys.index:
                another_mol = Chem.MolFromSmiles(dude_decoys.smiles.loc[decoy])
                fp_decoy = Chem.RDKFingerprint(another_mol)
                dude_actives_dude_decoys_tanimotos.append(DataStructs.TanimotoSimilarity(fp_active,fp_decoy)) 

            print("Done computing Tanomotos for DUD-E active", active, "and DUD-E decoy")

    csvs = glob("*_100_isomers.csv") # read in data for all 100 isomers of every active
    all_isomer_data = pd.DataFrame(index = range(0), columns = ['active', 'active_smi', 'active_MolWt', 'active_HBA', 'active_HBD', 'active_rotors', 'active_charge', 'active_TPSA', 'active_LabuteASA', 'isomer_smi', 'Tanimoto_similarity_to_active', 'isomer_MolWt', 'isomer_HBA', 'isomer_HBD', 'isomer_rotors', 'isomer_charge', 'isomer_TPSA', 'isomer_LabuteASA'])    
    for csv in csvs:
        temp = pd.read_csv(csv)
        all_isomer_data = pd.concat([all_isomer_data, temp], axis = 0)

    all_isomer_data.reset_index(inplace = True, drop = True)
    
###################################################### PLOT PROPERTIES FOR INDIVIDUAL TARGETS ########################################################

    _, p = mannwhitneyu(dude_actives_dude_decoys_tanimotos, all_isomer_data.Tanimoto_similarity_to_active)
    
###################################################### PLOT DUD-E ACTIVES - DUD-E DECOYS AND DUD-E ACTIVES - ISOMER DECOYS TANIMOTOS ########################################################

    plt.close()
    plt.title('Tanimoto similarities for ' + os.getcwd().split('/')[-1] + '\n' + 'p = ' + str(p), fontsize = 20)
    plt.hist(dude_actives_dude_decoys_tanimotos, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - experimental dude_decoys similarities')
    plt.hist(all_isomer_data.Tanimoto_similarity_to_active, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - isomerization dude_decoys similarities')
    plt.legend(fontsize=11)
    plt.xlabel('Tanimoto similarity', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.99, top=0.86, bottom=0.11)
    plt.savefig('Tanimoto_similarities.jpeg', dpi = 600)
    plt.close()
    
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS MOLECULAR WEIGHTS ########################################################

    plt.title('Molecular weight distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.MolWt, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.MolWt, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_MolWt, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('Molecular weight, [g/mol]', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_molecular_weights.jpeg', dpi = 600)
    plt.close()
    
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS HBAs ########################################################
    
    plt.title('HBA distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.HBA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.HBA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_HBA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('Number of HBA', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_hbas.jpeg', dpi = 600)
    plt.close()
    
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS HBDs ########################################################
    plt.title('HBD distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.HBD, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.HBD, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_HBD, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('Number of HBD', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_hbds.jpeg', dpi = 600)
    plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS ROTORS ########################################################

    plt.title('Rotatable bonds distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.rotors, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.rotors, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_rotors, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('Number of rotatable bonds', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_rotors.jpeg', dpi = 600)
    plt.close()
    
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS CHARGES ########################################################

    plt.title('Charge distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.charge, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.charge, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_charge, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('Charge, [e]', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_charges.jpeg', dpi = 600)
    plt.close()              
    
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS TPSAs ########################################################
        
    plt.title('TPSA distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.TPSA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.TPSA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_TPSA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('TPSA, [Å²]', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_tpsa.jpeg', dpi = 600)
    plt.close()        
        
###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS LABUTEASAs ########################################################

    plt.title('LabuteASA distributions for ' + os.getcwd().split('/')[-1], fontsize = 20)
    plt.hist(actives.LabuteASA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
    plt.hist(dude_decoys.LabuteASA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
    plt.hist(all_isomer_data.isomer_LabuteASA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
    plt.legend(fontsize=11)
    plt.xlabel('LabuteASA, [Å²]', fontsize = 18)
    plt.ylabel('Density', fontsize = 18)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
    plt.savefig('Tanimoto_LabuteASA.jpeg', dpi = 600)
    plt.close()                               

    os.chdir('..')

    return actives, dude_decoys, all_isomer_data, dude_actives_dude_decoys_tanimotos     
    
run = Parallel(n_jobs=20)(delayed(compare_ligs)(folders[i]) for i in range(0,len(folders))) # run processing for all 20 targets simultaneously and aggregate the data

actives = pd.concat([run[0][0],run[1][0],run[2][0],run[3][0],run[4][0],run[5][0],run[6][0],run[7][0],run[8][0],run[9][0],run[10][0],run[11][0],run[12][0],run[13][0],run[14][0],run[15][0],run[16][0],run[17][0],run[18][0],run[19][0]])

dude_decoys = pd.concat([run[0][1],run[1][1],run[2][1],run[3][1],run[4][1],run[5][1],run[6][1],run[7][1],run[8][1],run[9][1],run[10][1],run[11][1],run[12][1],run[13][1],run[14][1],run[15][1],run[16][1],run[17][1],run[18][1],run[19][1]])

all_isomer_data = pd.concat([run[0][2],run[1][2],run[2][2],run[3][2],run[4][2],run[5][2],run[6][2],run[7][2],run[8][2],run[9][2],run[10][2],run[11][2],run[12][2],run[13][2],run[14][2],run[15][2],run[16][2],run[17][2],run[18][2],run[19][2]])

dude_actives_dude_decoys_tanimotos = run[0][3] + run[1][3] + run[2][3] + run[3][3] + run[4][3] + run[5][3] + run[6][3] + run[7][3] + run[8][3] + run[9][3] + run[10][3] + run[11][3] + run[12][3] + run[13][3] + run[14][3] + run[15][3] + run[16][3] + run[17][3] + run[18][3] + run[19][3]


###################################################### PLOT PROPERTIES FOR ENTIRE DATA SET ########################################################

###################################################### PLOT DUD-E ACTIVES - DUD-E DECOYS AND DUD-E ACTIVES - ISOMER DECOYS TANIMOTOS ########################################################

_, p = mannwhitneyu(dude_actives_dude_decoys_tanimotos, all_isomer_data.Tanimoto_similarity_to_active)

plt.close()
plt.title('Tanimoto similarities' + '\n' + 'p = ' + str(p), fontsize = 20)
plt.hist(dude_actives_dude_decoys_tanimotos, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - experimental dude_decoys similarities')
plt.hist(all_isomer_data.Tanimoto_similarity_to_active, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives - isomerization dude_decoys similarities')
plt.legend(fontsize=11)
plt.xlabel('Tanimoto similarity', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.99, top=0.86, bottom=0.11)
plt.savefig('Tanimoto_similarities.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS MOLECULAR WEIGHTS ########################################################

plt.title('Molecular weight distributions', fontsize = 20)
plt.hist(actives.MolWt, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.MolWt, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_MolWt, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('Molecular weight, [g/mol]', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_molecular_weights.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS HBAs ########################################################

plt.title('HBA distributions', fontsize = 20)
plt.hist(actives.HBA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.HBA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_HBA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('Number of HBA', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_hbas.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS HBDs ########################################################

plt.title('HBD distributions for', fontsize = 20)
plt.hist(actives.HBD, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.HBD, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_HBD, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('Number of HBD', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_hbds.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS ROTORS ########################################################

plt.title('Rotatable bonds distributions', fontsize = 20)
plt.hist(actives.rotors, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.rotors, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_rotors, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('Number of rotatable bonds', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_rotors.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS CHARGES ########################################################

plt.title('Charge distributions', fontsize = 20)
plt.hist(actives.charge, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.charge, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_charge, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('Charge, [e]', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_charges.jpeg', dpi = 600)
plt.close()  

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS TPSAs ########################################################

plt.title('TPSA distributions', fontsize = 20)
plt.hist(actives.TPSA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.TPSA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_TPSA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('TPSA, [Å²]', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_tpsa.jpeg', dpi = 600)
plt.close()

###################################################### PLOT DUD-E ACTIVES, DUD-E DECOYS, AND ISOMER DECOYS LABUTEASAs ########################################################

plt.title('LabuteASA distributions', fontsize = 20)
plt.hist(actives.LabuteASA, density=True, color = 'blue', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental actives')
plt.hist(dude_decoys.LabuteASA, density=True, color = 'green', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'experimental dude_decoys')
plt.hist(all_isomer_data.isomer_LabuteASA, density=True, color = 'red', alpha = 0.5, histtype = 'step', linewidth = 3, label = 'isomerization dude_decoys')
plt.legend(fontsize=11)
plt.xlabel('LabuteASA, [Å²]', fontsize = 18)
plt.ylabel('Density', fontsize = 18)
plt.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.11)
plt.savefig('Tanimoto_LabuteASA.jpeg', dpi = 600)
plt.close()   

