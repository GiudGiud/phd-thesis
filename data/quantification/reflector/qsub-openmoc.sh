#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=null" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=infinite" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=degenerate" openmoc-reflector.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=null" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=infinite" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=degenerate" openmoc-reflector.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=null" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=infinite" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=degenerate" openmoc-reflector.pbs