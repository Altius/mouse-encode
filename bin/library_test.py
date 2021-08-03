#!/usr/bin/env python
import pandas as pd,numpy as np,os,sys
import files_from_lims,cell_distances

PEAK_DIR = '/net/seq/data/projects/ENCODE3_publications/mouse_dnase_recall/hotspot2-fdr0.05/peaks'
CC_DIR = '/net/seq/data/projects/ENCODE3_publications/mouse_dnase_recall/perBase/results'

def main(metadata,master_list_file,output_dir):
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    metadata = pd.read_table(metadata,header = 0)

    new_metadata = get_library_data(metadata)
    mfile = os.path.join(output_dir,'library-metadata.txt')
    new_metadata.to_csv(mfile,sep = '\t',index = True)

    cell_distances.main(mfile,master_list_file,output_dir)

def get_library_data(metadata):
    new_metdata = {'DS-ID':[],'Method':[],'Cut-Count-File':[],'Peak-File':[],
                   'Tissue':[],'Time':[],'Date':[],'Prep-Date':[]}
    for i,ds in enumerate(list(metadata['DS-ID'].values)):
        try:
            nmb = int(ds.lstrip('DS'))
        except ValueError:
            nmb = int(ds.lstrip('DS')[0:-1])
        library_data = files_from_lims.get_ds_libraries(nmb)
        sample_data = files_from_lims.get_sample_data(nmb)
        for sub_id,method,date in library_data:
            new_ds = sub_id
            print new_ds
            print library_data
            if new_ds in new_metdata['DS-ID']: continue
            cc_file = [os.path.join(CC_DIR,x)
                       for x in os.listdir(CC_DIR) if new_ds in x and x.endswith('starch')]
            pk_file = [os.path.join(PEAK_DIR,x)
                       for x in os.listdir(PEAK_DIR) if new_ds in x]
            if not cc_file or not pk_file: continue
            new_metdata['DS-ID'].append(new_ds)
            new_metdata['Date'].append(date)
            new_metdata['Prep-Date'].append(sample_data['prep_date'])
            new_metdata['Method'].append(method)
            new_metdata['Cut-Count-File'].append(cc_file[0])
            new_metdata['Peak-File'].append(pk_file[0])
            if i >= metadata.shape[0]:
                new_metdata['Tissue'].append('hindbrain')
                new_metdata['Time'].append('P0')
            else:
                new_metdata['Tissue'].append(metadata.iloc[i]['Tissue'])
                new_metdata['Time'].append(metadata.iloc[i]['Time'])

    new_metdata = pd.DataFrame.from_dict(new_metdata)
    new_metdata = new_metdata.set_index('DS-ID')
    return new_metdata

if __name__ == "__main__":
    main(*sys.argv[1:])