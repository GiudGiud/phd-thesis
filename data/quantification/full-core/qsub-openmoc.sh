#!/bin/bash

#qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=2 --clusterizer=null -i 10000" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=2 --clusterizer=infinite -i 10000" openmoc-full-core.pbs
#qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=2 --clusterizer=degenerate -i 10000" openmoc-full-core.pbs

#qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=8 --clusterizer=null -i 10000" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=8 --clusterizer=infinite -i 10000" openmoc-full-core.pbs
#qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=8 --clusterizer=degenerate -i 10000" openmoc-full-core.pbs

#qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=null -i 10000" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=infinite -i 10000" openmoc-full-core.pbs
#qsub -P moose  -v flags="-a 64 -s 0.05 -c 1E-7 --num-fine=70 --clusterizer=degenerate -i 10000" openmoc-full-core.pbs
