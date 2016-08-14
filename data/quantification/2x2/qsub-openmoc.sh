#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=null" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=infinite" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=degenerate" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=null" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=infinite" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=degenerate" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=null" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=infinite" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=degenerate" openmoc-2x2.pbs