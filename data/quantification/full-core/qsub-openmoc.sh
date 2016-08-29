#!/bin/bash

qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=null -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=infinite -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=2 --clusterizer=degenerate -i 200" openmoc-full-core.pbs

qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=null -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=infinite -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=8 --clusterizer=degenerate -i 200" openmoc-full-core.pbs

qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=40 --clusterizer=null -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=40 --clusterizer=infinite -i 200" openmoc-full-core.pbs
qsub -P moose -W depend=afterany:701890.service2 -v flags="-a 128 -s 0.05 --num-fine=40 --clusterizer=degenerate -i 200" openmoc-full-core.pbs
