#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=null" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=degenerate" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=lns" openmoc-2x2.pbs