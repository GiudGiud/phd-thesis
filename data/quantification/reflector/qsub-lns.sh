#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=lns" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=lns" openmoc-reflector.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --clusterizer=lns" openmoc-reflector.pbs