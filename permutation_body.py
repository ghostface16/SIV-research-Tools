from permutation_core import*
import pandas as pd
import numpy as np
from math import comb

def permutation_body(labels_unique, df_w_lab:pd.DataFrame, nsamp_array:dict, prop_to_permut:float, 
    average_size_minus_4:int, dim:int, size_loop:int, sep:str,
    labels_unique_indx, random_labels, pairwise_ks_matrix:np.array, 
    pKs_matrix:np.array, ana_loc:str, nproc:int, n_permut:int, pairwise:bool,
    count_ana_locs_1:int, partition_2:int, compute_random_Ks:bool):

    #print(df_w_lab)
    df_no_lab = df_w_lab.drop(df_w_lab.index[-1], axis=0)# df_w_lab.drop(index = 'Labels') 
    animal = df_w_lab.columns[0].split(sep)[0]
    nsamp_keys = list(nsamp_array.keys())
    labels_row = df_w_lab.iloc[-1,:]
    #print(labels_unique)

    #for ana_loc in labels_unique:
    ref_df = df_no_lab.loc[:,labels_row==ana_loc]
    ref_df_cols = ref_df.columns
    nsamp = nsamp_array[ana_loc]
    nsamp_minus_2 = nsamp - 2
    partition = int(comb(nsamp, 2) )
    pKs_array = np.zeros((n_permut,))
    n_seqs_to_perm = int(round(nsamp*prop_to_permut))  #####!!!!!######!!!!!!!***** 
    
    if not pairwise:
        w1 = float(nsamp_minus_2/(nsamp+average_size_minus_4))
        w2 = float(1-w1)
        rand_cols = None
        random_df = None
        pairwise_ks_matrix = None
        count_ana_locs = None
    else:
        count_ana_locs = dim - size_loop  #####!!!!!######!!!!!!!*****
        result = None

    for i in range(size_loop):
        if pairwise:
            ana_loc_to_compare = labels_unique[count_ana_locs]
            nsamp_2_key = nsamp_keys[count_ana_locs]
            nsamp_2 = nsamp_array[nsamp_2_key]
            partition_2 = int(comb(nsamp_2, 2))
            random_df = df_no_lab.loc[:,labels_row==ana_loc_to_compare]
            rand_cols = random_df.columns
            average_size_minus_4 = nsamp_2 - 4 
            #n_seqs_to_perm = round(nsamp_2*prop_to_permut)  #####!!!!!######!!!!!!!***** 
            w1 = float(nsamp_minus_2/(nsamp+average_size_minus_4))
            w2 = float(1-w1)
        
        result = permutation_core(n_permut=n_permut, pairwise=pairwise, ref_df=ref_df, labels_unique_indx=labels_unique_indx, 
                                  random_labels=random_labels, df_no_lab=df_no_lab, n_seqs_to_perm=n_seqs_to_perm, 
                                  ref_df_cols=ref_df_cols, nproc=nproc, count_ana_locs_1=count_ana_locs_1, w1=w1, w2=w2, 
                                  partition=partition, partition_2=partition_2, rand_cols=rand_cols, random_df=random_df, 
                                  pairwise_ks_matrix=pairwise_ks_matrix, count_ana_locs=count_ana_locs)
        if pairwise:
            count_ana_locs+=1

    while not result:
        pass

    if compute_random_Ks:
        #pKs_matrix = result[0]
        pKs_matrix = pd.DataFrame(result[0], columns=[ana_loc])
        return((animal, ana_loc, pKs_matrix))
    else:
        pKs_matrix = pd.DataFrame(result[0], columns=[ana_loc])
        return((animal, ana_loc, pKs_matrix))
        print('###WORK IT OUT')






