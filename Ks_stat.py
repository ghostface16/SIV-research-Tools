from math import floor, log, comb
import pandas as pd
import numpy as np

def Ks_stat(df:pd.DataFrame, pairwise:bool, grouping_index:int, size_Ki_samp=100, compute_random_Ks=True):

    df_cols = df.columns

    # Very application specific, must be replaced by code that takes a row as an input an sort sequences by location or other kinds of label
    # It should be performed by user and not implemented into code
    #labels_all = [x.split('.')[1] + '_' + x.split('.')[2] + '_' + x.split('.')[3] + '_' + x.split('.')[4] if len(x.split('.'))>5 else x.split('.')[1] + '_' + x.split('.')[2] + '_' + x.split('.')[3] if len(x.split('.'))>4 else x.split('.')[1] + '_' + x.split('.')[2] if len(x.split('.'))>3 else x.split('.')[1] if len(x.split('.'))>2 else x.split('.')[0] for x in df.columns]
    
    labels_all = df.iloc[grouping_index, :]

    labels_unique_and_count = np.unique(labels_all, return_counts=True)
    labels_count = labels_unique_and_count[1]
    labels_unique = labels_unique_and_count[0]

    #Should allow user to choose minimum number of seq per group for a group to be included in computations
    is_one = np.where(labels_count==1)[0]

    # even more application specific
    is_22617 = np.where(labels_unique=='>22617')[0]
    
    # This is fine we can keep it since we want to get rid of groups with less then X members
    if len(is_one)>0:
        labels_unique = np.delete(labels_unique, is_one)
        labels_count = np.delete(labels_count, is_one)

    # Even more application specific should be taken out too, users will have to deal with there own dataset specific issues
    if len(is_22617)>0:
        labels_unique = np.delete(labels_unique, is_22617 )
        labels_count = np.delete(labels_count, is_22617 )

    labels_df = pd.DataFrame({'Labels':labels_all}, index=df_cols).transpose()
    
    # since we're letting user figuring out grouping before running the test, this line isn't needed anymore
    #df_w_labels = pd.concat([df, labels_df], axis=0)
    df_w_labels = df
    print(df_w_labels)
    labels_row = df_w_labels.iloc[grouping_index,:]
    K_array = {}#np.zeros(len(labels_unique))
    nsamp_array = {}
    df_no_lab = df_w_labels.drop(df_w_labels.index[-1], axis=0)#index='Labels')
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
        Ks_df = pd.DataFrame(Ks_array, columns=labels_unique)
        return (K_array, nsamp_array, random_K_array, Ks_df, random_labels,
                df_w_labels,(labels_unique, labels_count), average_size)
    elif pairwise:
        return (pairwise_ks_matrix, K_array, nsamp_array, (labels_unique, labels_count), 
                df_w_labels)
    else:
        return (K_array, nsamp_array, (labels_unique, labels_count), df_w_labels)
