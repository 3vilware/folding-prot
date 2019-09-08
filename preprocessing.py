# This file is part of the OpenProtein project.
#
# @author Jeppe Hallgren
#
# For license information, please see the LICENSE file in the root directory.

import glob
import os.path
import os
import platform
import numpy as np
import h5py
from util import AA_ID_DICT, calculate_dihedral_angles, protein_id_to_str, get_structure_from_angles, \
    structure_to_backbone_atoms, write_to_pdb, calculate_dihedral_angles_over_minibatch, \
    get_backbone_positions_from_angular_prediction, encode_primary_string
import torch

MAX_SEQUENCE_LENGTH = 2000

def process_raw_data(use_gpu, force_pre_processing_overwrite=True):
    print("Starting pre-processing of raw data...")
    input_files = glob.glob("data/raw/*")
    print(input_files)
    input_files_filtered = filter_input_files(input_files)
    for file_path in input_files_filtered:
        if platform.system() is 'Windows':
            filename = file_path.split('\\')[-1]
        else:
            filename = file_path.split('/')[-1]
        preprocessed_file_name = "data/preprocessed/"+filename+".hdf5"

        # check if we should remove the any previously processed files
        if os.path.isfile(preprocessed_file_name):
            print("Preprocessed file for " + filename + " already exists.")
            if force_pre_processing_overwrite:
                print("force_pre_processing_overwrite flag set to True, overwriting old file...")
                os.remove(preprocessed_file_name)
            else:
                print("Skipping pre-processing for this file...")

        if not os.path.isfile(preprocessed_file_name):
            process_file(filename, preprocessed_file_name, use_gpu)
    print("Completed pre-processing.")

def read_protein_from_file(file_pointer):

        dict_ = {}
        _dssp_dict = {'L': 0, 'H': 1, 'B': 2, 'E': 3, 'G': 4, 'I': 5, 'T': 6, 'S': 7}
        _mask_dict = {'-': 0, '+': 1}

        while True:
            next_line = file_pointer.readline()
            print ("nextline:", next_line)
            if next_line == '[ID]\n':
                id_ = file_pointer.readline()[:-1]
                dict_.update({'id': id_})
            elif next_line == '[PRIMARY]\n':
                primary = encode_primary_string(file_pointer.readline()[:-1])
                dict_.update({'primary': primary})
            elif next_line == '[EVOLUTIONARY]\n':
                evolutionary = []
                for residue in range(21): evolutionary.append(
                    [float(step) for step in file_pointer.readline().split()])
                dict_.update({'evolutionary': evolutionary})
            elif next_line == '[SECONDARY]\n':
                secondary = list([_dssp_dict[dssp] for dssp in file_pointer.readline()[:-1]])
                dict_.update({'secondary': secondary})
            elif next_line == '[TERTIARY]\n':
                tertiary = []
                # 3 dimension
                for axis in range(3): 
                    tertiary.append([float(coord) for coord in file_pointer.readline().split()])
                dict_.update({'tertiary': tertiary})
            elif next_line == '[MASK]\n':
                mask = list([_mask_dict[aa] for aa in file_pointer.readline()[:-1]])
                dict_.update({'mask': mask})
            elif next_line == '\n':
                return dict_
            elif next_line == '':
                return None

def read_protein_from_file_ids(file_pointer):

        dict_ = {}
        _dssp_dict = {'L': 0, 'H': 1, 'B': 2, 'E': 3, 'G': 4, 'I': 5, 'T': 6, 'S': 7}
        _mask_dict = {'-': 0, '+': 1}

        while True:
            next_line = file_pointer.readline()
            print ("nextline:", next_line)
            if next_line == '[ID]\n':
                id_ = file_pointer.readline()[:-1]
                print("ID:", id_)
                id_ = id_.split('_')[0]
                f = open('data/protein_ids.txt', 'a+') # append mode
                f.write(id_ + '\n')

            elif next_line == '':
                f.close()
                return None


import time
def process_file(input_file, output_file, use_gpu):
    print("Processing raw data file", input_file)

    # create output file
    f = h5py.File(output_file, 'w')
    current_buffer_size = 1
    current_buffer_allocation = 0
    dset1 = f.create_dataset('primary',(current_buffer_size,MAX_SEQUENCE_LENGTH),maxshape=(None,MAX_SEQUENCE_LENGTH),dtype='int32')
    dset2 = f.create_dataset('tertiary',(current_buffer_size,MAX_SEQUENCE_LENGTH,9),maxshape=(None,MAX_SEQUENCE_LENGTH, 9),dtype='float')
    dset3 = f.create_dataset('mask',(current_buffer_size,MAX_SEQUENCE_LENGTH),maxshape=(None,MAX_SEQUENCE_LENGTH),dtype='uint8')

    input_file_pointer = open("data/raw/" + input_file, "r")

    while True:
        # while there's more proteins to process
        next_protein = read_protein_from_file(input_file_pointer)
        if next_protein is None:
            print("NEXT IS NONE")
            break
        print("ID PROT",next_protein["id"])
        sequence_length = len(next_protein['primary'])

        if sequence_length > MAX_SEQUENCE_LENGTH:
            print("Dropping protein as length too long:", sequence_length)
            continue

        if current_buffer_allocation >= current_buffer_size:
            current_buffer_size = current_buffer_size + 1
            dset1.resize((current_buffer_size,MAX_SEQUENCE_LENGTH))
            dset2.resize((current_buffer_size,MAX_SEQUENCE_LENGTH, 9))
            dset3.resize((current_buffer_size,MAX_SEQUENCE_LENGTH))

        primary_padded = np.zeros(MAX_SEQUENCE_LENGTH)
        tertiary_padded = np.zeros((9, MAX_SEQUENCE_LENGTH))
        mask_padded = np.zeros(MAX_SEQUENCE_LENGTH)

        # masking and padding here happens so that the stored dataset is of the same size. 
        # when the data is loaded in this padding is removed again. 

        # todo: here is the original tertiary structure (FIRST ACHIVEMENT)
        """
        # Uncomment
        primary_padded[:sequence_length] = next_protein['primary'] # padd with zeros max seq letters number 
        t_transposed = np.ravel(np.array(next_protein['tertiary']).T) # Flatten layer 
        t_reshaped = np.reshape(t_transposed, (sequence_length,9)).T
        """

        tertiary = np.array(next_protein['tertiary'])/100
        primary = np.array(next_protein['primary'])
        print("SHAPE TERTIARY", tertiary.shape)
        print("SHAPE PRIMARY", primary.shape)
        seq_aux_len = sequence_length
        cont = 0
        aux_dim = 3
        tertiary = (tertiary.T)
        print("LAST:", tertiary[-1], 0)
        print("LAST:", tertiary[-1], 1)
        print("LAST:", tertiary[-1], 2)

        while True:
            aux_dim = cont + 3
            print("\nCont", cont, "AUX", aux_dim)
            print(tertiary[cont:aux_dim, 0])
            #print(tertiary[cont:aux_dim, 1])
            #print(tertiary[cont:aux_dim, 2])
            cont += sequence_length
            x = input('continue...')

        print("PRIMARY SHAPE", np.shape(next_protein['primary']))

        """
        print("ORIGINAL PRIMARY",np.shape(next_protein['primary']), next_protein['primary'] ) # 
        print("MASKED PRIMARY",np.shape(next_protein['primary']), next_protein['primary'][32:] ) # 
        print("PRIMARY PADDED",np.shape(primary_padded), primary_padded ) # 
        print("ORIGINAL TERTIARY", next_protein['tertiary'] ) # np.shape is 3,726 like in document
        print("ORIGINAL TERTIARY SHAPE", np.shape(next_protein['tertiary'] ))

        print("[ Tertiary from doc transposed ]", np.shape(t_transposed))
        print("[ Tertiary from doc reshaped into 9 ]", np.shape(t_reshaped))
        # time.sleep(15000)
        #for i in t_reshaped:
        #    for e in i:
        #        print("Coord:", e)
        #time.sleep(15000)
        """


        tertiary_padded[:,:sequence_length] = t_reshaped
        mask_padded[:sequence_length] = next_protein['mask']

        print("MASK PADDED", np.shape(mask_padded), mask_padded)
        mask = torch.Tensor(mask_padded).type(dtype=torch.uint8)
        print("MASK TENSOR", np.shape(mask), mask)
        
        prim = torch.masked_select(torch.Tensor(primary_padded).type(dtype=torch.long), mask)
        print("PRIMARY SEQ INFO ON TENSOR CONVERSION:", np.shape(prim), ":\n", prim)
        pos = torch.masked_select(torch.Tensor(tertiary_padded), mask).view(9, -1).transpose(0, 1).unsqueeze(1) / 100 # wtf
        pos_aux = torch.masked_select(torch.Tensor(tertiary_padded), mask).view(9, -1)#.transpose(0, 1)#.unsqueeze(1) / 100 # wtf
        print("POS SHAPE:", np.shape(pos_aux))
        #print("POS after mask selected:", pos_aux)

        """
        if use_gpu:
            pos = pos.cuda()

        angles, batch_sizes = calculate_dihedral_angles_over_minibatch(pos, [len(prim)], use_gpu=use_gpu)

        tertiary, _ = get_backbone_positions_from_angular_prediction(angles, batch_sizes, use_gpu=use_gpu)
        tertiary = tertiary.squeeze(1)

        primary_padded = np.zeros(MAX_SEQUENCE_LENGTH)
        tertiary_padded = np.zeros((MAX_SEQUENCE_LENGTH, 9))


        length_after_mask_removed = len(prim)
        print("DEBUG:", prim)
        print("PRIM CPU.NUMPY():", prim)
        print("SHAPE TERTIARY (3, 768)
PRIMARY SHAPE (256,)
LENGTH AFTER MASK REMOVED:", length_after_mask_removed)
        primary_padded[:length_after_mask_removed] = prim.data.cpu().numpy()
        tertiary_padded[:length_after_mask_removed, :] = tertiary.data.cpu().numpy()
        mask_padded = np.zeros(MAX_SEQUENCE_LENGTH)
        mask_padded[:length_after_mask_removed] = np.ones(length_after_mask_removed)
        print("MASK FULLY PADDED", np.shape(mask_padded), mask_padded)

        dset1[current_buffer_allocation] = primary_padded
        dset2[current_buffer_allocation] = tertiary_padded
        dset3[current_buffer_allocation] = mask_padded
        current_buffer_allocation += 1
        
        print("Buffer allocation", current_buffer_allocation)
        print("\nDSET final primary ",np.shape(primary_padded),":", primary_padded)
        print("\nDSET final tertiary",np.shape(tertiary_padded),":", tertiary_padded)
        print("Final Structure")
        # x = input('Enter para continuar...')
        """
        break # to run the first protein
    print("Wrote output to", current_buffer_allocation, "proteins to", output_file)


def filter_input_files(input_files):
    disallowed_file_endings = (".gitignore", ".DS_Store")
    return list(filter(lambda x: not x.endswith(disallowed_file_endings), input_files))


#process_file('sample.txt', 'output/sample.h5', False)
#process_file('1vbk.txt', 'output/sample.h5', False)
#process_file('training_30', 'output/sample.h5', False)

