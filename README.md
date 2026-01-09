<p align="center">
<img width="6990" height="5200" alt="image" src="https://github.com/user-attachments/assets/1cd0829c-4438-4294-aa30-3c875674a3c2" />
</p>


This page contains the code used to generate and analyse isomers for 1640 active molecules for 20 targets from the DUD-E set along with the original DUD-E data for those targets (actives_final.ism and decoys_final.ism). The original DUD-E publication can be found at https://pubs.acs.org/doi/10.1021/jm300687e. 

**The work is described in detail** [here](https://chemrxiv.org/engage/chemrxiv/article-details/695b7900083c11e4a1deea81). **Also see** [this preprint](https://chemrxiv.org/engage/chemrxiv/article-details/68737a05728bf9025e3ae9c2).

To run it, make sure to install [Open Babel 3.1.1](https://openbabel.org/) and [MAYGEN](https://link.springer.com/article/10.1186/s13321-021-00529-9); Open Babel should be in your PATH.

First, cd into the input_data_and_code folder and generate isomers for every active in the set.


```
cd /path/to/input_data_and_code # adjust this properly 
for i in cxcr4 comt fabp4 pur2 glcm sahh pygm hs90a hxk4 ada mcr nram pa2ga fak1 hivint nos1 rock1 grik1 mapk2 def
do
cd $i
python ../generate_isomers.py
cd ..
done
```

Then, run isomer searches. In the code, you can provide SMARTS patterns for molecular features to avoid, and specify acceptable margins for the number of hydrogen bond acceptors (HBA), donors (HBD), and rotors. 


```
python compare_ligs.py
```

Finally, compute DUD-E active - DUD-E active, DUD-E decoy - DUD-E decoy, and isomer - isomer Tanimoto similarities. Note the last step can take a lot of RAM. All code herein is provided free of charge and without any warranty as an example; modify according to your needs and resources.


```
python compare_actives_dude_decoys_and_isomer_decoys.py
```
