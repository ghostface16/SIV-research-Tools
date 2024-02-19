from Ks_multiprocessing import*

#give full path to csv
df_file = sys.argv[2]
df = pd.read_csv(df_file, index_col=0)
nrep_par = int(sys.argv[1])

result = run_Kstar_perm(df, pairwise=False, compute_random_Ks=True, 
                        size_Ki_samp=1000,  n_permut=1000, prop_to_permut=0.2, nrep_par=nrep_par)