#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=null" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=degenerate" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=lns" openmoc-assm-3.1-20BPs.pbs