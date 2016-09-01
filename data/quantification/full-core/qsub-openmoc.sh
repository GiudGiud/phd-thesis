#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=null -i 50" openmoc-full-core.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=infinite -i 50" openmoc-full-core.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=degenerate -i 50" openmoc-full-core.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=null -i 50" openmoc-full-core.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=infinite -i 50" openmoc-full-core.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=degenerate -i 50" openmoc-full-core.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=null -i 50" openmoc-full-core.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=infinite -i 50" openmoc-full-core.pbs
qsub -P moose  -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=degenerate -i 50" openmoc-full-core.pbs
