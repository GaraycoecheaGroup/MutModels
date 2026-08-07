# this will run if dirichlet_multinomial is run as a module 
# so if the package is executed as a tool
# can put argparsing etc. here


def background_correct(data,background):
    """
    subtract background from counts, then remove nonnegative and round to counts
    """
    data = data.sub(background,axis=0)
    data[data<0] = 0
    data = data.round().astype(int)
    return data

def select_substitutions(data,substitutions):
    selection = []
    for sub in substitutions:
        cur = [i for i in data.index.values if sub in  i]

        selection+= (cur)

    # select the relevant substitutions without changing index order of selected 
    data = data[data.index.isin(selection)]
    return data

if __name__ == "__main__":
    # for now I run it directly to test it out
    from mutmodels.common.parsers import read_matrix,align
    import pandas as pd

    from mutmodels.dirichlet_multinomial.api import dm_two_condition_SBS
    # figure 4 test data
    # matrix_fn = "./test_data/fig4.SBS96.all"
    # bg_fn = './test_data/WT_UNT_background.SBS96.tsv'

    # figure 2 test data
    matrix_fn = "./test_data/fig2.SBS96.all"
    bg_fn = './test_data/fig2.UNT_background.SBS96.tsv'

    # figure 4 samples
    # g1 = ['JG0159n','JG0159o','JG0159p','JG0159q','JG0159r','JG0159s','JG0159t']
    # g2 = ['JG0182b','JG0182c','JG0182d','JG0182e']

    # figure 2 failed to fit samples:
    fdur_sg2 = ['JG0137b','JG0137c','JG0137d','JG0137e']
    fdur_sg3 = ['JG0138b','JG0138c','JG0138d','JG0138e']

    # figure 2 samples
    g2 = ['JG0139b','JG0139c','JG0139d','JG0139e']
    g1 = ['JG0140b','JG0140c','JG0140d','JG0140e']

    # substitutions to include
    substitutions = ['T>G','T>C','T>A']

    # test the api function
    res = dm_two_condition_SBS(matrix_fn,g1,g2,bg_fn=bg_fn,substitutions=substitutions,rng=1234)

    print('api test:')
    print(res['p_value'])
    print(res['TVD'])
    print()

    # test of the thing that failed to fit
    # res = dm_two_condition_SBS(matrix_fn,fdur_sg2,fdur_sg3,bg_fn = bg_fn,substitutions=substitutions,rng=1234)
    # print('failtofit test:')


    # print(res['p_value'])
    # print(res['TVD'])
    # print(res['a0_A'])
    # print(res['a0_B'])