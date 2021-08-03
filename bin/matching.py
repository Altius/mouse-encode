#!/usr/bin/env python
import numpy as np,pandas as pd
import subprocess

SEED_NUMBER = 2000
SEED = np.random.RandomState(SEED_NUMBER)

PSEUDOCOUNT = 1
MATCHING_TOLERANCE = 0.05

def main(bed_file,master_file,matrix_file,tss_file,output_file,tolerance = 0.005,per_sample = 1):
    bed_subset = get_subset(bed_file,master_file)
    master_values = get_signal(matrix_file)
    not_subset = np.setdiff1d(np.arange(len(master_values)),bed_subset)

    case_distances = get_tss_distances(bed_file,tss_file)
    control_distances = get_tss_distances(master_file,tss_file)[not_subset]


    case_values = master_values[bed_subset]
    control_values = master_values[not_subset]

    distance_options = quantile_match(case_distances,control_distances,tolerance)
    value_options = quantile_match(case_values,control_values,tolerance)

    background_set = []
    for i in xrange(len(case_values)):
        opts = np.intersect1d(distance_options[i],value_options[i])
        background_set.extend(list(not_subset[SEED.choice(opts,replace = False,size = per_sample)]))

    write_subset(background_set,master_file,output_file)



def quantile_match(values,vector,percentage):
    arg_order = np.argsort(vector)
    offset = int(percentage * len(vector))
    rlocations = np.searchsorted(vector[arg_order], values, side='right')
    llocations = np.searchsorted(vector[arg_order], values, side='left')
    locations = (rlocations/2. + llocations/2.).astype(int)
    blocations = np.minimum(llocations,np.maximum(0,locations - offset))
    tlocations = np.maximum(rlocations,np.minimum(len(vector),locations + offset))
    return [arg_order[blocations[i]:tlocations[i]] for i in xrange(len(values))]


def write_subset(background_index,master_file,output_file):
    count = 0
    subset = set([x for x in background_index])
    fout = open(output_file,'w')
    with open(master_file) as f1:
        for line in f1:
            if count in subset: fout.write(line)
            count += 1
    fout.close()

def get_tss_distances(bed_file,tss_file):
    p1 = subprocess.Popen(['closest-features','--dist','--closest',bed_file,tss_file],stdout = subprocess.PIPE)
    values = []
    for line in iter(p1.stdout.readline,''):
        bd1,bd2,x  = line.split('|')
        values.append(abs(int(float(x))))
    p1.stdout.close()
    return np.array(values)

def get_signal(matrix_file):
    dt = pd.read_table(matrix_file,index_col = 0,header = 0)
    return np.mean(dt.values,axis = 1)

def get_subset(bed1,master_file):
    p1 = subprocess.Popen(['bedmap','--indicator',bed1,master_file],stdout = subprocess.PIPE)
    vls = []
    for line in iter(p1.stdout.readline,''):
        vls.append(int(float(line.strip())))
    p1.stdout.close()
    return np.where(np.array(vls))[0]


def nearest_neighbor_matching(background_values, sample_values, tolerance=MATCHING_TOLERANCE, number_per_sample=1,
                              log_transform=False, pseudocount=1):
    if log_transform:
        background_values = np.log(background_values + 1)
        sample_values = np.log(sample_values + 1)

    sorted_ind_order = np.argsort(background_values)
    if number_per_sample < 1:
        sub_index = SEED.choice(len(background_values), int(number_per_sample * len(background_values)), replace=False)
        sorted_ind_order = sorted_ind_order[sub_index]



    #### if there are not enough peaks within the tolerance, choose the closest ones
    bottom_locations = np.maximum(0, np.minimum(bottom_locations, locations - number_per_sample / 2))
    top_locations = np.minimum(len(background_values) - 1, np.maximum(top_locations, locations + number_per_sample / 2))

    sample = set([])
    for i, l_mid in enumerate(locations):
        l_high = top_locations[i]
        l_low = bottom_locations[i]
        if l_high == l_low and l_mid != len(sorted_ind_order):
            sample.add(sorted_ind_order[l_mid])
            ### this would be odd if reached
            continue
        potential = [x for x in xrange(l_low, l_high) if sorted_ind_order[x] not in sample]
        if len(potential) == 0: continue
        indices = SEED.choice(potential, min(number_per_sample, len(potential)))

        sample.update([sorted_ind_order[r] for r in indices])
    return np.array(list(sample), dtype=int)