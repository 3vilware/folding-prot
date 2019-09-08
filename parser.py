import requests
AA_ID_DICT = {'ALA': 1, 'CYS': 2, 'ASP': 3, 'GLU': 4, 'PHE': 5, 'GLY': 6, 'HIS': 7, 'ILE': 8, 'LYS': 9,
              'LEU': 10, 'MET': 11, 'ASN': 12, 'PRO': 13, 'GLN': 14, 'ARG': 15, 'SER': 16, 'THR': 17,
              'VAL': 18, 'TRP': 19,'TYR': 20}

AA_ID_DICT_SINGLE = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
                    'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
                    'V': 18, 'W': 19,'Y': 20}

visited = {}
coords = []

def get_fasta(pdb_id):
    begin = False
    full_sequence = ""
    result = requests.get(
            'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList={}&compressionType=uncompressed'.format(pdb_id))
    for line in str(result.text).splitlines():
        if line.startswith('>') and not begin:
            begin = True
        elif line.startswith('>') and begin:
            break
        else:
            full_sequence = full_sequence + line
    print(list(full_sequence))
    return list(full_sequence)

def get_pdb(pdb_id):
    result = requests.get(
            'https://files.rcsb.org/view/{}.pdb'.format(pdb_id))
    #print(result.text)
    return result.text

def encode_sequence(seq):
    aa_encoded = []        
    for aa in seq:
        aa = aa.replace(' ','')
        for k, v in AA_ID_DICT_SINGLE.items():
            k = k.replace(' ','')
            if k == aa:
                aa_encoded.append(v)
    print(aa_encoded)
    return aa_encoded

def read_pdb(pdb_content):
    last_aa = ""
    aa_list = []
    for line in pdb_content.splitlines():
        if line[0:4] == 'ATOM':
            if last_aa == line[17:21]:
                pass
            else:
                #print("CAMBIANDO AMINIACIDO")
                last_aa = line[17:21]
                aa_list.append(line[17:21])
            #print (line[17:21], "X = ", line[32:39], "Y = ", line[40:48], "Z = ", line[49:56])
    return aa_list

def struct_pdb_to_tensor(pdb_content, sequence): # string sequence (not encoded)
    """
    leer cada atomo y agregarlo a un array diferente por cada x, y ,z para cada aa
    ejemplo:
        MET: M
        np.array([ [ |[x1,x2,x3..],[y1,y2,y3...], [z1,z2,z3...] ]| ]) para cada aminoacido en la cadena 
    """
    last_aa = ""
    aa_list = []
    for line in pdb_content.splitlines():
        if line[0:4] == 'ATOM':
            if last_aa == line[17:21]:
                pass
            else:
                #print("CAMBIANDO AMINIACIDO")
                last_aa = line[17:21]
                aa_list.append(line[17:21])
            #print (line[17:21], "X = ", line[32:39], "Y = ", line[40:48], "Z = ", line[49:56])
    return aa_list


def read_pdb_from_file(pdb_path):
    last_aa = ""
    aa_list = []
    for line in open(pdb_path):
        if line[0:4] == 'ATOM':
            if last_aa == line[17:21]:
                pass
            else:
                #print("CAMBIANDO AMINIACIDO")
                last_aa = line[17:21]
                aa_list.append(line[17:21])
            # print (line[17:21], "X = ", line[32:39], "Y = ", line[40:48], "Z = ", line[49:56])
    return aa_list



aa_seq = get_fasta('1VBK')
print(len(aa_seq))
encoded_aa = encode_sequence(aa_seq)
print(len(encoded_aa))
read_pdb(get_pdb('1VBK'))
