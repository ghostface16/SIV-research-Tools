from Ks_multiprocessing import*
import pandas as pd
import sys

# Give full path to csv
df_file = sys.argv[2]
df = pd.read_csv(df_file, index_col=0)
# Number of cores to use
nrep_par = int(sys.argv[1])

if len(sys.argv)>3:
    # Where to dump result CSVs, default is working dir
    csv_path = sys.argv[3]

result = run_Kstar_perm(df, nrep_par=nrep_par, pairwise=False, compute_random_Ks=True, csv_dump=True, 
                        grouping_index= -1, size_Ki_samp=1000,  n_permut=1000, prop_to_permut=0.2, sep="_")
