#!/usr/bin/env python
import sys,files_from_lims,pandas as pd

def main(metadata_file,output_file):

    dt = pd.read_table(metadata_file,index_col = 0,header = 0)
    for vl in dt.index.values:
        agv = dt.loc[vl,'AG']
        if type(agv) is str: agv = int(agv.lstrip('AG'))
        print agv
        cc_file,pk_file = files_from_lims.get_file_paths(agv)
        dt.loc[vl,'Cut-Count-File'] = cc_file
        dt.loc[vl,'Peak-File'] = pk_file

    if 'DS_plus' in dt.columns.values:
        dt['DS-ID'] = dt['DS_plus'].values
        dt = dt.set_index(keys = 'DS_plus',drop = False)
    dt.to_csv(output_file,sep = '\t')


if __name__ == "__main__":
    main(*sys.argv[1:])