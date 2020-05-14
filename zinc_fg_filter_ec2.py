import rdkit
from rdkit import Chem
import os, glob
from time import time
from datetime import datetime

mlc_dict = {}

mlc_dict['ketones'] = '[#6&!$(C#N)][CX3](=[OX1])[#6&!$(C#N)]'
mlc_dict['primary_amines'] = '[NX3H2][c,CX4&!$([CX4]([NH2])[O,N,S,P])]'
mlc_dict['secondary_amines'] = '[NX3H1]([c,CX4&!$([CX4]([NX3H1])[O,S,N,P])])[c,CX4&!$([CX4]([NX3H1])[O,S,N,P])]'
mlc_dict['tertiary_amines'] = '[NX3]([c,CX4&!$([CX4]([NX3])[O,S,N,P])])([c,CX4&!$([CX4]([NX3])[O,S,N,P])])[c,CX4&!$([CX4]([NX3])[O,S,N,P])]'
mlc_dict['alkyl_halides'] = '[CX4][Cl,Br,I]'
mlc_dict['aryl_halides'] = 'a[Cl,Br,I]'
mlc_dict['organozincs'] = '[Zn][#6]'
mlc_dict['nitriles'] = '[#1,#6&!$([CX3]=[OX1,SX1])][CX2]#[NX1]'
mlc_dict['acyl_chlorides'] = '[#1,#6][CX3](=[OX1])Cl'
mlc_dict['carbamates'] = '[#6][OX2][CX3](=[OX1])[NH2,NH1,NX3]'
mlc_dict['sulfonic_acid_esters'] = '[#6][Sv6X4](=[OX1])(=[OX1])[OX2][#6]'
mlc_dict['boronic_acid_derivs'] = '[Bv3X3]([#6])([!#6])[!#6]'
mlc_dict['alkynes'] = '[CX2]#[CX2]'
mlc_dict['arenes'] = 'c1ccccc1'
mlc_dict['benzyl_chlorides'] = 'a[CX4]Cl'
mlc_dict['aromatic_5m_monoheterocycles'] = '[a!#6]1cccc1'
mlc_dict['aromatic_5m_diheterocycles'] = '[$([a!#6]1[a!#6]ccc1),$([a!#6]1c[a!#6]cc1)]'
mlc_dict['aromatic_5m_triheterocycles'] = '[$([a!#6]1[a!#6][a!#6]cc1),$([a!#6]1[a!#6]c[a!#6]c1)]'
mlc_dict['aromatic_6m_monoheterocycles'] = '[a!c]1ccccc1'
mlc_dict['aromatic_6m_diheterocycles'] = '[$([a!c]1[a!c]cccc1),$([a!c]1c[a!c]ccc1),$([a!c]1cc[a!c]cc1)]'
mlc_dict['aromatic_6m_triheterocycles'] = '[$([a!c]1[a!c][a!c]ccc1),$([a!c]1[a!c]c[a!c]cc1),$([a!c]1c[a!c]c[a!c]c1)]'
mlc_dict['alkenes'] = '[!#8!#7!#16][CX3]([!#8!#7!#16])=[CX3]([!#8!#7!#16])[!#8!#7!#16]'
mlc_dict['alkenyl_halides'] = '[Cl,Br,I][CX3]=[CX3]'
mlc_dict['hydroxylamines'] = '[N&!$([NH1]([CX3]=[OX1,SX1,NX2])[OH1])&!$(N([OH1])([CX3]=[OX1,SX1,NX2])([CX3]=[OX1,SX1,NX2]))]([#1,#6])([#1,#6])[OH1]'
mlc_dict['all_C_silanes'] = '[Si](-[#6])(-[#6])(-[#6])-[#6]'

curr_path = os.getcwd()
smi_folder = 'zinc_remaining'
fg_dict = mlc_dict
path = os.path.join(curr_path, smi_folder)
files = glob.glob(f'{path}/*')
ts = time()
start_time = datetime.now()
print(f'Running FG filter. Started at {start_time}.')
for file in files:
    head, tail = os.path.split(file)
    zinc_set = tail[-5]
    dir_name = f'{zinc_set}_mlc_fgs'
    new_dir_name = os.path.join(curr_path, dir_name)
    os.makedirs(new_dir_name)
    results = []
    mol_tups = []
    with open(file) as infile:
        print(f'Processing file: {tail}')
        for line in infile:
            line = line.split()
            if not len(line) == 0:
                smi = line[0]
                if not Chem.MolFromSmiles(smi) == None:
                    mol_tups.append((smi, Chem.MolFromSmiles(smi)))
        print(f'Mol writing for set {zinc_set} complete.')
        print(len(mol_tups), 'smiles read')
    for name, smarts in fg_dict.items():
        if not Chem.MolFromSmarts(smarts) == None:
            smol = Chem.MolFromSmarts(smarts)
            match_count = 0
            match_list = []
            for tup in mol_tups:
                smiles = tup[0]
                mol = tup[1]
                try:
                    mol.HasSubstructMatch(smol)
                except:
                    continue
                if mol.HasSubstructMatch(smol) == True:
                    match_count +=1
                    match_list.append(smiles)
            with open(f'{new_dir_name}/{zinc_set}_{name}.txt', 'w') as outfile:       
                for entry in match_list:
                    outfile.write(entry + '\n')
            results.append(f'{name}: {match_count}')
        else:
            continue
    results = sorted(results)
    with open(f'{new_dir_name}/{zinc_set}_report.txt', 'w') as report_file:
        for entry in results:
            report_file.write(entry + '\n')
    print(f'Filtering set {zinc_set} complete.')
te = time()
ttm = round(((te-ts)/60), 3)
end_time = datetime.now()
print('Processing complete.')
print(f'Finished at {end_time}.')
print(f'Full processing time: {ttm} min.')