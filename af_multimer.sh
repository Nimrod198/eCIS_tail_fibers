#!/bin/bash

#load GPU software stuff
module load cuda/12.4.1

#single chain

#sbatch \
#--gres gpu:a100-1-10 \
#--time=23:0:0 \
#--mem=24G \
#--wrap="colabfold_batch --num-recycle 10  --recycle-early-stop-tolerance=0.5 --stop-at-score 90 --model-type alphafold2_ptm --rank ptmscore ./make_fastas/$1 out_$1"

#complex

#sbatch \
#--gres gpu:a100-1-10 \
#--time=36:0:0 \
#--mem=24G \
#--wrap="colabfold_batch --num-recycle 20 --recycle-early-stop-tolerance=0.5 --stop-at-score 90 --model-type alphafold2_multimer_v3 --rank multimer ./make_fastas/$1 out_$1"

#complex in a30 gpu

#--partition=puffin \
#--gres gpu:a30:2 \


sbatch \
--partition=puffin \
--gres gpu:a30:3 \
--time=120:0:0 \
--cpus-per-task=2 \
--mem-per-cpu=25G \
--wrap="/sci/labs/asafle/alexlevylab/icore-data/tools/local_collabfold/v1.5.0/localcolabfold/colabfold-conda/bin/colabfold_batch --num-recycle 5 --num-models 1 --recycle-early-stop-tolerance=0.5 --sort-queries-by length --model-type alphafold2_multimer_v3 --rank multimer $1 $2"

