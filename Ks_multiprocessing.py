from math import floor, log
import pandas as pd
import numpy as np
from permuted_Ks_multiprocessing import*
import os

       
################################# pKs permutation statistic  ##############################################################################

def run_Kstar_perm(df:pd.DataFrame,  nrep_par:int, pairwise:bool, compute_random_Ks:bool, csv_dump:bool,
                        size_Ki_samp=100,  n_permut=1000, prop_to_permut=0.2, csv_path=None):

    results = permuted_Ks_multiprocessing(df=df, pairwise=pairwise, compute_random_Ks=compute_random_Ks, 
                        size_Ki_samp=size_Ki_samp,  n_permut=n_permut, prop_to_permut=prop_to_permut, 
                        nrep_par=nrep_par)

    if csv_dump:
        if not csv_path:
            csv_path = os.getcwd()

        ident = results[-1]
        ks_name = ident + '_Ks.csv'
        ks_path = os.path.join(csv_path, ks_name)
        pks_name = ident + '_pKs.csv'
        pks_path = os.path.join(csv_path, pks_name)
        Ks_df = results[0]
        pKs_df = results[2]
        pKs_df.to_csv(f'{pks_path}')
        Ks_df.to_csv(f'{ks_path}')

    return(pKs_df, Ks_df)

##########################################################################################################


def compute_p_val(Ks_stat:pd.DataFrame, pKs_matrix:pd.DataFrame, signif_thresh=0.001, fdr_alpha = 0.1):
    
    H = 0
    p_val_dict = {}
    alpha = signif_thresh/2
    one_minus_alpha = 1 - alpha
    trials = pKs_matrix.shape[0]
    upper_limit_ks = float(Ks_stat.quantile(one_minus_alpha)[0])
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
