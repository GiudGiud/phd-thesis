#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=null" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=degenerate" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=lns" openmoc-refl.pbs