from Kaplan_distance import*
import pandas as pd
import numpy as np
#from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor as Pool

def permutation_core(n_permut:int, pairwise:bool, ref_df:pd.DataFrame, labels_unique_indx, random_labels,
                    df_no_lab:pd.DataFrame, n_seqs_to_perm:int, ref_df_cols, nproc:int,
                     count_ana_locs_1, w1:float, w2:float, partition, partition_2):

    if not pairwise:
        pKs_array = np.zeros((n_permut,))

    for permutation in range(n_permut):
        #picking one of the random populations generated by Ks_stat()  
        if not pairwise:    
            random_samp = np.random.choice(labels_unique_indx,1)
            rand_cols = random_labels[f'{int(random_samp)}']
            random_df = df_no_lab[rand_cols]

        #picking ref and random sequences to permute
        random_perm_col_index = np.random.choice(rand_cols, n_seqs_to_perm, replace=False)
        ref_perm_col_index = np.random.choice(ref_df_cols, n_seqs_to_perm, replace=False)
        #permutting ref and random sequences
        permutted_random_df = random_df.drop(random_perm_col_index, axis=1)
        permutted_random_df = pd.concat([permutted_random_df, ref_df[ref_perm_col_index]],axis=1)
        permutted_ref_df = ref_df.drop(ref_perm_col_index, axis=1)
        permutted_ref_df = pd.concat([permutted_ref_df, random_df[random_perm_col_index]],axis=1)
        #permutted random and ref columns
        permutted_ref_df_cols = permutted_ref_df.columns
        permutted_random_df_cols = permutted_random_df.columns
        permutation_core_args = [permutted_ref_df, permutted_random_df]

        if nproc>1:
            with Pool(max_workers=nproc) as pool_2:
                #task = pool_2.starmap_async(Kaplan_distance, permutation_core_args)
                #results = pool.map(Kaplan_distance, permutation_core_arg)
                results_0 = [pool_2.submit(pKaplan_distance, arg) for arg in permutation_core_args]
                #task.close()    #Close Pool and let all the processes complete
                #task.join()
                #results = task.result() 
                results = [result.result() for result in results_0]
        else:
          results = [pKaplan_distance(df) for df in permutation_core_args]
        
        previous = results[0][0]
        previous_2 = results[1][0]
        pK = previous/partition
        rpK = previous_2/partition_2
        pKw = pK*w1
        rpKw = rpK*w2
        
        if pairwise:
            pairwise_ks_matrix[count_ana_locs_1-1, count_ana_locs, permutation] = pKw + rpKw
        else:
            pKs_array[permutation] = pKw + rpKw

    if pairwise:
        return((pairwise_ks_matrix[count_ana_locs_1-1, count_ana_locs,:],))
    else:
        #pKs_matrix =
        pKs_matrix = pd.DataFrame(pKs_array, columns=[ana_loc])
        return ((pKs_matrix[ana_loc],))


