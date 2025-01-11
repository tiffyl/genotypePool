#!/usr/bin/env python3

import pandas as pd
import altair as alt
import sys

alt.data_transformers.disable_max_rows()

## INPUT
gt_inputlist = sys.argv[1].replace(' ','').strip('[]').split(',')

## FUNCTIONS
def idCellLine(pool_result_gtfile):
    genotypes = pd.read_table(pool_result_gtfile).query('`agree.crick` >= 5 & `agree.watson` >= 5').rename(columns=(lambda x: x.replace('.', "_")))
    genotypes['frac_agree_crick'] = genotypes['agree_crick'] / genotypes['all_pos_crick']
    genotypes['frac_agree_watson'] = genotypes['agree_watson'] / genotypes['all_pos_watson']

    #Each line in the file is a unique pair of cellline to single cell
    maxid = pd.DataFrame({'crick': genotypes.groupby(by=['cell_ID'])['frac_agree_crick'].idxmax(), 'watson':genotypes.groupby(by=['cell_ID'])['frac_agree_watson'].idxmax()}).query('crick == watson')
    cells = genotypes[genotypes.index.isin(maxid['crick'])].query('frac_agree_crick >= 0.8 | frac_agree_watson >= 0.8')
    
    chr=genotypes['chr'].unique()[0]
    return(cells[['cell_ID', 'sample', 'frac_agree_crick', 'frac_agree_watson']].rename(columns=(lambda x: x if x in ["chr", "cell_ID"] else chr+"_"+x)))
    
## SCRIPT
cells_df = pd.DataFrame()

for gtfile in gt_inputlist:
    cells = idCellLine(gtfile)
    
    if cells_df.empty:
        cells_df = cells
    else:
        cells_df = cells_df.merge(cells, how="outer", on='cell_ID')

cells_df = cells_df.set_index('cell_ID').filter(like='sample').reset_index()
cellline_df = cells_df.melt(id_vars='cell_ID').drop(columns=['variable']).dropna().rename(columns={'value':'cellline'}).groupby(['cell_ID', 'cellline']).size().reset_index(name='count').query('count > 1')

cells_df.to_csv('allChrsCellLineAssign.tsv', sep="\t", index=False)
cellline_df[['cell_ID', 'cellline']].to_csv('scCellLineAssign.tsv', sep="\t", index=False, header=False)
pd.DataFrame(cellline_df.value_counts('cellline')).reset_index().to_csv('countCellLineAssign.tsv', sep="\t", index=False, header=False)