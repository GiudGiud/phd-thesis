#!/bin/bash

qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=2 --clusterizer=combined --num-clusters=30 -i 35" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=8 --clusterizer=combined --num-clusters=30 -i 35" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --clusterizer=combined --num-clusters=30 -i 35" openmoc-full-core.pbs