from pathlib import Path
import os
import glob
import pandas as pd
from pandas import DataFrame
import numpy as np

def collect_energies(work_dir: str) -> DataFrame:
    work_dir = Path(work_dir).absolute()
    print(f"Extract energies from mdout files in {work_dir}")
    files = glob.glob(str(work_dir / 'min_fr*.mdout'))
    #print(files)
    order = [f.split('/')[-1] for f in files]
    order = [int(f.replace('.mdout', '').split('_')[-1]) for f in order]
    order = np.argsort(order)
    files = np.array(files)[order]
    results = []

    for file in files:
        name = file.split('/')[-1].replace('.mdout', '')
        f = open(file, 'r')
        lines = f.readlines()
        #print(len(lines))
        #print(file)
        header_line_id = [i for i in range(
            len(lines)) if 'FINAL RESULTS' in lines[i]][0]
        tot_E_id = header_line_id + 5
        bond_angle_dih_id = header_line_id + 7
        vdw_eel_gb_id = header_line_id + 8
        one_four_id = header_line_id + 9

        totE = lines[tot_E_id].split()[1]           
        bond = lines[bond_angle_dih_id].split()[2]  
        angle = lines[bond_angle_dih_id].split()[5] 
        dihed = lines[bond_angle_dih_id].split()[8] 
        vdw = lines[vdw_eel_gb_id].split()[2]       
        eel = lines[vdw_eel_gb_id].split()[5]       
        egb = lines[vdw_eel_gb_id].split()[8]       
        one_four_vdw = lines[one_four_id].split()[3]
        one_four_eel = lines[one_four_id].split()[7]
        
        results.append(
            [name, totE, bond, angle,
             dihed, vdw, eel, egb,
             one_four_vdw, one_four_eel
             ]
        )
        f.close()
    results = pd.DataFrame(results, columns=[
        'name', 'ENERGY', 'BOND', 'ANGLE',
        'DIHED', 'VDW', 'EEL', 'EGB',
        '1-4 VDW', '1-4 EEL'
    ])
    print("All energies in kcal/mol")
    return results


# 1. Collect energies of minimized frames
pdb = 'aB13' # output name
energies = collect_energies(work_dir='mdout/')
energies.to_csv(pdb + '_energies.csv', index=False)
