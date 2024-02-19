from math import floor, log
import pandas as pd
import numpy as np
from permuted_Ks_multiprocessing import*


def Ks_stat(df:pd.DataFrame, pairwise:bool, size_Ki_samp=100, compute_random_Ks=True):

    df_cols = df.columns
    labels_all = [x.split('.')[1] + '_' + x.split('.')[2] + '_' + x.split('.')[3] + '_' + x.split('.')[4] if len(x.split('.'))>5 else x.split('.')[1] + '_' + x.split('.')[2] + '_' + x.split('.')[3] if len(x.split('.'))>4 else x.split('.')[1] + '_' + x.split('.')[2] if len(x.split('.'))>3 else x.split('.')[1] if len(x.split('.'))>2 else x.split('.')[0] for x in df.columns]
    labels_unique_and_count = np.unique(labels_all, return_counts=True)
    labels_count = labels_unique_and_count[1]
    labels_unique = labels_unique_and_count[0]
    is_one = np.where(labels_count==1)[0]
    is_22617 = np.where(labels_unique=='>22617')[0]
    
    if len(is_one)>0:
        labels_unique = np.delete(labels_unique, is_one)
        labels_count = np.delete(labels_count, is_one)
    if len(is_22617)>0:
        labels_unique = np.delete(labels_unique, is_22617 )
        labels_count = np.delete(labels_count, is_22617 )

    labels_df = pd.DataFrame({'Labels':labels_all}, index=df_cols).transpose()
    df_w_labels = pd.concat([df, labels_df], axis=0)
    labels_row = df_w_labels.iloc[-1,:]
    K_array = {}#np.zeros(len(labels_unique))
    nsamp_array = {}
    df_no_lab = df_w_labels.drop(index='Labels')
    random_labels = {}
    
    if compute_random_Ks:
        Ks_array = {}
        average_size = int(labels_count.mean())
        random_K_array = np.zeros((size_Ki_samp,))
        random_partition = comb(average_size, 2)
        average_size_minus_4 = average_size - 4
        
        for iters in range(size_Ki_samp):
            random_cols = np.random.choice(df_cols, average_size, replace=False)
            random_labels[f'{iters}'] = random_cols
            random_ref = df_no_lab[random_cols]
            count, previous = 0, 0
            random_ref_cols = random_ref.columns
            
            for seq_1 in random_ref_cols:
                for seq_2 in random_ref_cols[count:]:
                    comparison = random_ref[seq_1]==random_ref[seq_2]
                    D = log(comparison.sum()+1)
                    previous = previous + D
                
                count+=1
        
            random_K = previous/random_partition 
            random_K_array[iters] = random_K
        
    count_2 = 0
    
    if pairwise:  
        dim = len(labels_count)
        pairwise_ks_matrix = np.zeros((dim, dim))#####!!!!!######!!!!!!!***** 
    else:
        size_loop = 1
 
    for ana_loc, nsamp in zip(labels_unique, labels_count):
        ref_df = df_no_lab.loc[:,labels_row==ana_loc]
        ref_df_cols = ref_df.columns
        count, previous = 0, 0
        partition = comb(nsamp, 2)
                
        for seq_1 in ref_df_cols:
            for seq_2 in ref_df_cols[count:]:
                comparison = ref_df[seq_1]==ref_df[seq_2]
                D = log(comparison.sum()+1)
                previous = previous + D

            count+=1
        
        K = previous/partition 
        nsamp_array[ana_loc] = nsamp

        if compute_random_Ks:
            K_array[ana_loc] = K
            w1 = (nsamp-2)/(nsamp+average_size_minus_4)
            w2 = 1-w1
            Kiw = K*w1
            rKw = random_K_array*w2
            Ks_array[ana_loc] = rKw + Kiw
        else:
            K_array[ana_loc] = K

        count_2+=1
        
    if pairwise:
        pairwise_count = 0
        k_list = list(K_array.values())
        
        for val_k, nsamp in zip(k_list, labels_count): 
            count_pairwise_2 = pairwise_count +1
            nsamp_minus_4 = nsamp - 4
            
            for val_k_2, nsamp_2 in zip(k_list[count_pairwise_2:], labels_count[count_pairwise_2:]):
                w1 = (nsamp-2)/(nsamp_minus_4+nsamp_2)
                w2 = 1-w1
                Ks = (val_k*w1) + (val_k_2*w2)
                pairwise_ks_matrix[pairwise_count, count_pairwise_2] = Ks
                count_pairwise_2+=1
            
            pairwise_count+=1
    
    if compute_random_Ks:
        return (K_array, nsamp_array, random_K_array, Ks_array, random_labels,
                df_w_labels,(labels_unique, labels_count), average_size)
    elif pairwise:
        return (pairwise_ks_matrix, K_array, nsamp_array, (labels_unique, labels_count), 
                df_w_labels)
    else:
        return (K_array, nsamp_array, (labels_unique, labels_count), df_w_labels)


################################# pKs permutation statistic  ##############################################################################

def run_Kstar_perm(df:pd.DataFrame, pairwise:bool, compute_random_Ks:bool, 
                        size_Ki_samp=100,  n_permut=1000, prop_to_permut=0.2)

    results = permuted_Ks_multiprocessing(df=df, pairwise=pairwise, compute_random_Ks=compute_random_Ks, 
                        size_Ki_samp=size_Ki_samp,  n_permut=n_permut, prop_to_permut=prop_to_permut)
    return(results[0])

##########################################################################################################


def compute_p_val(Ks_stat:pd.DataFrame, pKs_matrix:pd.DataFrame, signif_thresh=0.001, fdr_alpha = 0.1):
    
    H = 0
    p_val_dict = {}
    alpha = signif_thresh/2
    one_minus_alpha = 1 - alpha
    trials = pKs_matrix.shape[0]
    upper_limit_ks = float(Ks_stat.quantile(one_minus_alpha))
    print(upper_limit_ks)
    
    for pks in pKs_matrix.columns:
        ana_loc = pKs_matrix[pks] 
        #upper_limit_ks = Ks_stat[pks].quantile(one_minus_alpha)
        lower_limit_pks = ana_loc.quantile(alpha)
        to_check_1 = lower_limit_pks <= ana_loc 
        to_use_pks = ana_loc[to_check_1]
        to_check = to_use_pks <= upper_limit_ks
        H = to_check.sum()
        p_val = H/trials
        p_val_dict[pks] = p_val
    
    all_p_vals = list(p_val_dict.values())
    print(all_p_vals)
    adj_p_vals = fdrcorrection(all_p_vals, alpha=fdr_alpha, method='indep', is_sorted=False)[1]
    p_val_dict_adj = {key:adj_p_val for key,adj_p_val in zip(p_val_dict.keys(),adj_p_vals)}
    
    return p_val_dict_adj#p_val_dict

################################################################################################

def compute_pairwise_p_val(Ks_stat:np.array, pKs_matrix:np.array, labels_dict:dict):
    
    H = 0
    p_val_dict = {}
    trials = np.shape(pKs_matrix)[2]
    lines, columns = np.shape(Ks_stat)[0], np.shape(Ks_stat)[0]
    labels_list = labels_dict
    
    for line in range(lines):
        line_list = []
        
        for column in range(columns):
            if column>line:
                ana_loc = pKs_matrix[line,column,:] 
                to_check = ana_loc <= Ks_stat[line,column] 
                H = to_check.sum()
                p_val = H/trials
            else:
                p_val = 0

            line_list.append(p_val)
            
        line_name = labels_list[line] 
        p_val_dict[line_name] = line_list
        
    #fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False)
        
    return p_val_dict