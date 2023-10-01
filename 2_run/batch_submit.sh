#!/bin/bash

###################################################################
# Author : Matt Bumstead, McMaster University, 2016               #
# Script to prepare a batch submission on:  wobbie.sharcnet.ca    #
#
#  - run this from dir:   morphologies/2_run/
#
#  Execute using:
#
#   ./batch_submit.sh > log/output_for_morphologies_.dat &
#
###################################################################

## ------ Date
#date="year-month-day"
date="2017-04-01"

## ------ Estimated Time
#TM="7d"
TM="1h"

## ------ Estimated RAM needed
#mppRAM="1G"
mppRAM="1G"

## ------ Quantifiers for sqsub -
shape="Circle"
size_range="36"
sample="hardbox"

#jname="polygon_numberofpolygons_boundary"
jname=$shape"_"$size_range"n_"$sample


outpath="/scratch/bumstema/Data/PolygonConfigurations/morphologies/source/3_data/data/jname/meta/simulation_log"


##------ ------ ------
for run in {36..36}
do
#echo sqsub --lazy -r $TM -q serial --mpp=1G -j $shape"_"$size_range"n_"$sample"_"$run -o log/output.$run.dat ./RunMorphologies $run &>/dev/null
sqsub --lazy -r $TM -q serial --mpp=$mppRAM -j $jname"_"$run -o ../3_data/data/jname/meta/simulation_log/log_output_$run.dat ./RunMorphologies $run &> /dev/null
sleep 0.1
done

#echo sqsub --lazy -r $TM -q serial --mpp=1G -j $shape"_"$size_range"n_"$sample"_"$run -o log/output.$run.dat ./RunMorphologies $run &>/dev/null
sqjobs
