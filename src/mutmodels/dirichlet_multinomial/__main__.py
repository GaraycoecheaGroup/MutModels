# this will run if dirichlet_multinomial is run as a module 
# so if the package is executed as a tool
# can put argparsing etc. here

if __name__ == "__main__":
    
    # for now I run it directly to test it out 

    from mutmodels.common.parsers import read_matrix,align
    import pandas as pd

    # get example matrix
    matrix_fn = "./test_data/fig4.SBS96.all"
    matrix = read_matrix(matrix_fn)

    # get example background signal
    bg_fn = './test_data/WT_UNT_background.SBS96.tsv'
    bg = pd.read_csv(bg_fn,sep='\t',index_col=0).squeeze("columns")
    
    print(bg)

    print(matrix.head())

    print(matrix.head(60)).index.values

