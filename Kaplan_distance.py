from math import floor, log
import pandas as pd
import numpy as np

def pKaplan_distance(compartment_df:pd.DataFrame):

    count, previous = 0, 0
    compartment_df_cols = compartment_df.columns

    for seq_1 in compartment_df_cols:
        to_use_1 = compartment_df[seq_1]
        len_to_use_1 = len(to_use_1.shape)

        if len_to_use_1>1:
            to_use_1 = to_use_1.iloc[:,0]

        for seq_2 in compartment_df_cols[count:]:
            to_use_2 = compartment_df[seq_2]
            len_to_use_2 = len(to_use_2.shape)

            if len_to_use_2>1:
                to_use_2 = to_use_2.iloc[:,0]

            comparison = to_use_1==to_use_2    
            D = log(comparison.sum()+1)
            previous = previous + D

        count+=1 

    return((previous,))