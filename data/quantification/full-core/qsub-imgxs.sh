#!/bin/bash

qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=2" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=4" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=8" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=16" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=32" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=pinch --num-clusters=64" openmoc-full-core.pbs

qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=2" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=4" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=8" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=16" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=32" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=combined --num-clusters=64" openmoc-full-core.pbs
