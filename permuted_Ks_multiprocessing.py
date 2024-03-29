############################ write code for permutation statistic pKs ##############################################################################

#ONLY ncpu >= # ana_loc for the moment

####################################################################################################################################################

#from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor as Pool
import pandas as pd
import numpy as np
from permutation_body import*
from Ks_stat import*
#import sys

def permuted_Ks_multiprocessing(df:pd.DataFrame, nrep_par:int, pairwise:bool, compute_random_Ks:bool, sep:str,
                                grouping_index:int, size_Ki_samp=100,  n_permut=1000, prop_to_permut=0.1):
    
    Ks_stat_results = Ks_stat(df=df, size_Ki_samp=size_Ki_samp, compute_random_Ks=compute_random_Ks, pairwise=pairwise, 
                                grouping_index=grouping_index)
    labels_unique = Ks_stat_results[-2][0]
    labels_unique_count = Ks_stat_results[-2][1]

    if pairwise:  
        dim = len(labels_unique_count)
        pairwise_ks_matrix = np.zeros((dim, dim, n_permut))#####!!!!!######!!!!!!!***** 
        size_loop = dim - 1
        count_ana_locs_1 = 1 #####!!!!!######!!!!!!!***** 
        df_w_lab = Ks_stat_results[-1]
        nsamp_array = Ks_stat_results[2]
        average_size_minus_4 = None
        labels_unique_indx = None
        random_labels = None
        pKs_matrix = None
        partition_2 = None
    else:
        random_labels = Ks_stat_results[-4]
        average_size = Ks_stat_results[-1]
        partition_2 = comb(average_size, 2)
        size_loop = 1
        pKs_matrix = {}
        average_size_minus_4 = average_size - 4
        labels_unique_indx = np.array(list(range(size_Ki_samp)))
        df_w_lab = Ks_stat_results[-3]
        nsamp_array = Ks_stat_results[1]
        dim = None
        pairwise_ks_matrix = None
        count_ana_locs_1 = None
        #n_seqs_to_perm = round(average_size*prop_to_permut)#####!!!!!######!!!!!!!***** 
    
    # Introducing multiprocessing
    n_rep_par = nrep_par
    n_ALs = len(labels_unique)

    with Pool(max_workers=n_rep_par) as pool:
        permutation_body_args = []

        if n_rep_par==n_ALs:
            nproc = 1   
            lucky_unique_count = None 
            
        elif n_rep_par>n_ALs:  
            ratio = n_rep_par/n_ALs

            if n_rep_par%n_ALs==0:
                nproc = int(ratio*n_ALs)
                lucky_unique_count = None
            else:
                floored_ratio = floor(ratio)
                #left = n_ALs-(floored_ratio*n_ALs)
                left = n_rep_par-floor((ratio-floored_ratio)*n_rep_par)
                prob = [1/n_ALs]*n_ALs
                lucky = np.random.choice(a=labels_unique, p=prob, size=left, replace=True)
                lucky_unique_count = np.unique(lucky, return_counts=True)
                nproc = int(floored_ratio)

        for ana_loc in labels_unique:
            if lucky_unique_count:
                if ana_loc in lucky_unique_count[0]:
                    where_lucky = np.where(lucky_unique_count[0]==ana_loc)
                    nproc = int(nproc + lucky_unique_count[1][where_lucky])

            args = (labels_unique, df_w_lab, nsamp_array, 
                    prop_to_permut, average_size_minus_4, dim, size_loop, sep,
                    labels_unique_indx, random_labels, pairwise_ks_matrix, 
                    pKs_matrix, ana_loc, nproc, n_permut, pairwise, count_ana_locs_1, 
                    partition_2, compute_random_Ks)      
            
            permutation_body_args.append(args)  
            
            if pairwise:
                count_ana_locs_1+=1
                size_loop-=1

        final_results_1 = [pool.submit(permutation_body, *arg) for arg in permutation_body_args]
        final_results = [result.result()[2] for result in final_results_1]

        
     #STILL NEED TO WORKOUT THE PAIRWISE result       
    if pairwise:
        pairwise_ks_matrix = pd.concat(final_results, axis=1)
        ident = [result.result()[0] for result in final_results_1][0]
        return(Ks_stat_results[0], pairwise_ks_matrix, Ks_stat_results[-2], ident)
    if compute_random_Ks:
        #result_list = {i[1]:i[2] for i in final_results}
        pKs_matrix = pd.concat(final_results, axis=1)
        #pd.DataFrame(result_dict)
        ident = [result.result()[0] for result in final_results_1][0]
        labels_n_counts = Ks_stat_results[-2]
        #labels = labels_n_counts[0]
        Ks_df = Ks_stat_results[3]#pd.DataFrame(Ks_stat_results[3], columns=labels)
        return(Ks_df, labels_n_counts, pKs_matrix, ident)

    else:        
        return(Ks_stat_results[0], pKs_matrix)



