#!/bin/bash --login


#$ -pe smp.pe 5

module load tools/gcc/cmake/3.23.0
module load libs/gcc/nlopt/2.6.2
module load apps/gcc/R/4.3.1

Rscript "/mnt/iusers01/ja01/m13477cc/scratch/cvLASSO.R"