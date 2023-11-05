from tqdm import tqdm
from glob import glob
import os 
def fsave(list_results, file, addn=False):
    if addn == True:
        list_results = [i+'\n' for i in list_results]
    with open(file, 'w') as f:
        f.writelines(list_results)
        f.close()

def fread(file, removen=True, remove1strow=False):
    with open(file, 'r') as f:
        contents = f.readlines()
        f.close()
    if removen == True:
        contents = [i.replace('\n','') for i in contents]
    if remove1strow == True:
        contents.pop(0)
    return contents
all_index = [] # to store all index from AFDB and PDB
# UniprotID, PDB_file_name, start_resi, end_resi, resolution
# add AFDB index to all_index
AFDB_resolution = 3.5 # treat AFDB as 3.5A resolution 
all_domain_list = os.listdir('AlphaFoldDB-domain-INDEX/')
for domain in tqdm(all_domain_list):
    domain_structure_prefix = domain.rsplit('-',1)[0] # in windows system
    uID = domain.split('-')[1]
    # print(domain_structure_prefix)
    contents = fread(f'AlphaFoldDB-domain-INDEX/{domain}')
    for line in contents:
        if ',' not in line:
            domain_NO, start_res, end_res = line.split(':')[0], line.split(':')[1].split('-')[0], line.split('-')[-1] # This step count the end as the structural 
            all_index.append(f'{uID},{domain_structure_prefix}-{domain_NO},{start_res},{end_res},{AFDB_resolution}')
        else:
            line_sub_domains = line.split(':')[1].split(',')
            domain_NO = line.split(':')[0]
            count_ = 0
            for sub_domain in line_sub_domains:
                start_res, end_res = sub_domain.split('-')[0], sub_domain.split('-')[-1]
                all_index.append(f'{uID},{domain_structure_prefix}-{domain_NO}-{count_},{start_res},{end_res},{AFDB_resolution}')
                count_ += 1
# add PDB index to all_index
Human_Uniprot_table = fread('PDB-INDEX/uniprot-organism__Homo+sapiens+(Human)+[9606]_.tab', remove1strow=True)
Human_UniprotIDs = [i.split()[0] for i in Human_Uniprot_table]
PDB_chain_UID_index = fread('PDB-INDEX/pdb_chain_uniprot.csv', remove1strow=True)
PDB_resolution_table = fread('PDB-INDEX/total_resolution_index.csv', remove1strow=True)

for uID in tqdm(Human_UniprotIDs):
    for PDB_chain in PDB_chain_UID_index:
        if f',{uID},' in PDB_chain:
            PDBid, chain, start_res, end_res = PDB_chain.split(',')[0],PDB_chain.split(',')[1],PDB_chain.split(',')[-2],PDB_chain.split(',')[-1]
            for resolution_line in PDB_resolution_table:
                if PDBid.upper()+',' in resolution_line:
                    if 'NMR,' not in resolution_line and 'SOLUTION SCATTERING' not in resolution_line:
                        resolution = resolution_line.split(',')[-2]
                    else:
                        resolution = 2.5 # 3A for NMR
                    all_index.append(f'{uID},{PDBid}_{chain},{start_res},{end_res},{resolution}')
                    # print(f'{uID},{PDBid}_{chain},{start_res},{end_res},{resolution}')

fsave(all_index, 'all_index3.5-nmr2.5_sub_domain.csv', addn=True)
