#!/bin/bash

qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=null" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=degenerate" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=lns" openmoc-full-core.pbs