#!/bin/bash


set -o pipefail

master_file=$1 #/home/jvierstra/proj/ftd/results.slurm/regions.bed
motif_file=$2
conservation_file=$3
pk_file=$4
node_number=$5
tss_file=$6
output_dir=$7

echo $@

code_dir=$(realpath $(dirname $0))

####
rm -r ${output_dir}/conservation

mkdir -p ${output_dir}
mkdir -p ${output_dir}/tf_beds
mkdir -p ${output_dir}/conservation
mkdir -p ${output_dir}/logs



cat <<__SCRIPT__ > ${output_dir}/slurm.write_tf_bed
#!/bin/bash
#
#SBATCH --output=${output_dir}/%J.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

set -e -o pipefail

$code_dir/write_loop_beds.py $master_file $tss_file $pk_file $node_number $output_dir/tf_beds
__SCRIPT__


JOB0=$(sbatch --job-name=write_tf.bed \
        ${output_dir}/slurm.write_tf_bed)
echo $JOB0


### next create and array job...
cat <<__SCRIPT__ > ${output_dir}/slurm.run_conservation
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
set -e -o pipefail

INPUT_FILE=\`ls ${output_dir}/tf_beds | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`

echo \$INPUT_FILE
$code_dir/get_motif_conservation.py ${output_dir}/tf_beds/\$INPUT_FILE \$INPUT_FILE $motif_file $conservation_file > \
$output_dir/conservation/\$INPUT_FILE.cons
__SCRIPT__

cat <<__SCRIPT__ > ${output_dir}/slurm.combine
#!/usr/bin/env python
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=1G
#SBATCH --cpus-per-task=8
import os,sys
dir = '${output_dir}/conservation'
weights = 0
bases = 0
for fl in os.listdir(dir):
    c,b = open(os.path.join(dir,fl)).readlines()[0].strip().split()
    cons = float(c)
    bp = int(b)
    weights += cons*bp
    bases += bp
sys.stdout.write('%f\n' %(weights/float(bases)))

__SCRIPT__


cat <<__SCRIPT__ > ${output_dir}/slurm.make_run_conservation
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=1G
#SBATCH --cpus-per-task=8

JOB1=\$(sbatch --export=ALL --job-name=run_conservation \
--array=1-\$(ls $output_dir/tf_beds | wc -l ) \
${output_dir}/slurm.run_conservation)

echo \$JOB1

JOB2=\$(sbatch --export=ALL --depend=afterok:\${JOB1##* } \
--job-name=run_conservation -o $output_dir/average_conservation.txt \
${output_dir}/slurm.combine)

echo \$JOB2

__SCRIPT__

JOB1=$(sbatch --export=ALL --depend=afterok:${JOB0##* } \
	--job-name=run_conservation \
	${output_dir}/slurm.make_run_conservation)


echo $JOB1
### then wait on that to parse