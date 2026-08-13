# this will run if dirichlet_multinomial is run as a module 
# so if the package is executed as a tool
# can put argparsing etc. here

if __name__ == "__main__":
    # # for now I run it directly to test it out

    from mutmodels.common.parsers import read_matrix,align
    import pandas as pd

    from mutmodels import dm_twosample_SBS,dm_onesample_SBS
    from mutmodels.dirichlet_multinomial.api import select_substitutions
    import mutmodels.dirichlet_multinomial.core as dm
    
    #############################
    # one-condition testcase(s) #
    #############################

    # read the matrix
    # matrix_fn = '.'
    cosmic_fn = '/home/newuser/Documents/data/reference_data/cosmic/COSMIC_catalogue-signatures_SBS96_v3.5/COSMIC_v3.5_SBS_GRCh38.txt'

    cosmic = pd.read_csv(cosmic_fn,sep='\t',index_col=0)
    cosmic.index.name = 'MutationType'

    # background corrected counts
    corr_matrix_fn = './test_data/fig3.corrected.SBS96.all'
    corr_matrix = read_matrix(corr_matrix_fn)
    # print(corr_matrix)

    # # read the reference signature
    ref_sig_fn = './test_data/cosmic_v3.5_GRCh38_SBS30.tsv'
    ref_sig = align(pd.read_csv(ref_sig_fn,sep='\t',index_col=0).squeeze("columns"))
    # print(ref_sig)

    # select substitutions, renormalise reference proportions
    ref_sig_CtoT = select_substitutions(ref_sig,['C>T'])
    ref_sig_CtoT = ref_sig_CtoT/ref_sig_CtoT.sum() # renormalized

    count_data = select_substitutions(corr_matrix,['C>T'])

    # print(ref_sig_CtoT.sum())

    # get counts
    obs_counts = count_data.values.T
    # print(obs_counts.shape)

    # perform the test
    res = dm.dm_onesample(obs_counts,ref_sig_CtoT.values,studentize=False,rng=1234,n_bootstraps=1000)
    print(res['p_value'])
    print(res['TVD'])

    new_res = dm_onesample_SBS(corr_matrix_fn,ref_sig_fn,substitutions=['C>T'],matrix_type='SBS96',n_bootstraps=1000,
                               studentize=False,rng=1234)
    print('now using api:')
    print(new_res['p_value'])
    print(new_res['TVD'])

    #############################
    # two-condition testcase(s) #
    #############################

    # from mutmodels.dirichlet_multinomial.api import dm_two_condition_SBS
    # figure 4 test data
    # matrix_fn = "./test_data/fig4.SBS96.all"
    # bg_fn = './test_data/fig4.WT_UNT_background.SBS96.tsv'

    # figure 2 test data
    matrix_fn = "./test_data/fig2.SBS96.all"
    bg_fn = './test_data/fig2.UNT_background.SBS96.tsv'

    # # figure 4 samples
    # g1 = ['JG0159n','JG0159o','JG0159p','JG0159q','JG0159r','JG0159s','JG0159t']
    # g2 = ['JG0182b','JG0182c','JG0182d','JG0182e']

    # # figure 2 comparing two fdur concentrations:
    # g1 = ['JG0138b','JG0138c','JG0138d','JG0138e']
    # g2 = ['JG0139b','JG0139c','JG0139d','JG0139e']

    # # figure 2 samples
    g2 = ['JG0139b','JG0139c','JG0139d','JG0139e']
    g1 = ['JG0140b','JG0140c','JG0140d','JG0140e']

    # # substitutions to include
    substitutions = ['T>G','T>C','T>A']

    # # test the api function
    res = dm_twosample_SBS(matrix_fn,g1,g2,bg_fn=bg_fn,substitutions=substitutions,rng=1234,)

    print('two condition api test:')
    print(res['p_value'])
    print(res['TVD'])
    print()

