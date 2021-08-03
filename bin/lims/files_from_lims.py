#!/usr/bin/env python
import requests,json,sys
from stamlims_api.rest import setup_api
from stamlims_api import lims

api = setup_api()

API_URL = 'https://lims.altiusinstitute.org/api'
LIMS_API_TOKEN='78c65586bb31d571dfe4523d1e5773984c4269da'
API_HEADERS = {'Authorization': 'Token {}'.format(LIMS_API_TOKEN)}

# TODO API scraping will change once LIMS endpoint is created. Please excuse the mess.

def get_tag_aggregations(tag_id,pg_size = 100):
    aggregations = []
    ds_ids = []
    js_data = api.get_single_result(url_addition='aggregation/?tag_id=%d&page_size=%d' %(tag_id,pg_size))
    js_data = json.loads(json.dumps(js_data))
    aggregations = [x['id'] for x in js_data['results']]
    ds_ids = ['DS'+str(x['sample_number'])+x['library_sublibrary'] for x in js_data['results']]
    while js_data['next'] is not None:
        js_data = api.get_single_result(url = js_data['next'])
        js_data = json.loads(json.dumps(js_data))
        ds_ids.extend(['DS'+str(x['sample_number'])+x['library_sublibrary'] for x in js_data['results']])
        aggregations.extend([x['id'] for x in js_data['results']])
    return aggregations,ds_ids

def get_ds_aggregations(ds_id,pg_size = 100):
    aggregations = []
    ds_ids = []
    js_data = api.get_single_result(url_addition='aggregation/?tag_id=%d&page_size=%d' %(tag_id,pg_size))
    js_data = json.loads(json.dumps(js_data))
    aggregations = [x['id'] for x in js_data['results']]
    ds_ids = ['DS'+str(x['sample_number'])+x['library_sublibrary'] for x in js_data['results']]
    while js_data['next'] is not None:
        js_data = api.get_single_result(url = js_data['next'])
        js_data = json.loads(json.dumps(js_data))
        ds_ids.extend(['DS'+str(x['sample_number'])+x['library_sublibrary'] for x in js_data['results']])
        aggregations.extend([x['id'] for x in js_data['results']])
    return aggregations,ds_ids

def get_file_paths(agg_id):
    # use aggregation id to get paths to cutcounts and peaks
    get_peaks = "{}/file/?format=json&object_content_type=125&object_id={}&purpose__slug=hotspot-peaks".format(
        API_URL, agg_id)
    get_counts = "{}/file/?object_id={}&object_content_type=125&purpose__slug=cutcounts-starch".format(
        API_URL, agg_id)
    try:
        v = requests.get(get_counts, headers=API_HEADERS).json()
        cutcounts = v['results'][0]['path']
        p = requests.get(get_peaks, headers=API_HEADERS).json()
        peaks = p['results'][0]['path']
        return cutcounts, peaks
    except Exception as e:
        print  "Error found: {}. Make sure peak and cutcount files exist for AG{}.".format(e, agg_id)
        return None, None

def main(tag_id,output_file):
    aggregations,ds_ids = get_tag_aggregations(tag_id)
    cut_counts_files = []
    peaks = []
    for ag in aggregations:
        cc,pk = get_file_paths(ag)
        cut_counts_files.append(cc)
        peaks.append(pk)
    lines = ['\t'.join(['DS-ID','Aggregation-ID','Cut-Count-File','Peak-File'])]
    for i,a in enumerate(aggregations):
        if cut_counts_files[i] is None: continue
        lines.append('\t'.join([str(x) for x in [ds_ids[i],aggregations[i],cut_counts_files[i],
                                                peaks[i]]]))
    with open(output_file,'w') as f1: f1.write('\n'.join(lines)+'\n')

if __name__ == "__main__":
    ag_id = int(sys.argv[1])
    output_file = sys.argv[2]
    main(ag_id,output_file)