#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=lns" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=lns" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=lns" openmoc-assm-1.6.pbs