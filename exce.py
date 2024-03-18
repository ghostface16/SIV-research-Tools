from Ks_multiprocessing import*
import pandas as pd
import sys

#give full path to csv
df_file = sys.argv[2]
df = pd.read_csv(df_file, index_col=0)
nrep_par = int(sys.argv[1])

if len(sys.argv)>3:
    csv_path = sys.argv[3]

result = run_Kstar_perm(df, nrep_par=nrep_par, pairwise=False, compute_random_Ks=True, csv_dump=True, g
                        rouping_index= -2, size_Ki_samp=1000,  n_permut=1000, prop_to_permut=0.2)
