#!/bin/bash

qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --num-clusters=2" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --num-clusters=4" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --num-clusters=6" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --num-clusters=8" openmoc-full-core.pbs
qsub -P moose -v flags="-a 64 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --num-clusters=10" openmoc-full-core.pbs