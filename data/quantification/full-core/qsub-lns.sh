#!/bin/bash

qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=2 --clusterizer=lns -i 10000" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=8 --clusterizer=lns -i 10000" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=lns -i 10000" openmoc-full-core.pbs