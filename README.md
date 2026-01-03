<img width="671" height="499" alt="image" src="https://github.com/user-attachments/assets/3129c9c6-77aa-4554-bac4-2ee32d0bf407" />





This page contains the code used to generate and analyse isomers for 1640 active molecules for 20 targets from the DUD-E set along with the original DUD-E data for those targets (actives_final.ism decoys_final.ism). The original DUD-E publication can be found at https://pubs.acs.org/doi/10.1021/jm300687e. 

To run it, make sure to properly install [Open Babel 3.1.1](https://openbabel.org/) and [MAYGEN](https://link.springer.com/article/10.1186/s13321-021-00529-9). 

First, generate isomers for every active in the set.

```
for i in cxcr4 comt fabp4 pur2 glcm sahh pygm hs90a hxk4 ada mcr nram pa2ga fak1 hivint nos1 rock1 grik1 mapk2 def
do
cd $i
python ../generate_isomers.py
cd ..
done
```
