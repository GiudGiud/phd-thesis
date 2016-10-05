#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=agglomerative --transformer=hmm --num-clusters=2" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=agglomerative --transformer=hmm --num-clusters=4" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=agglomerative --transformer=hmm --num-clusters=6" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=agglomerative --transformer=hmm --num-clusters=8" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=agglomerative --transformer=hmm --num-clusters=10" openmoc-assm-1.6.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=kmeans --transformer=hmm --num-clusters=2" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=kmeans --transformer=hmm --num-clusters=4" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=kmeans --transformer=hmm --num-clusters=6" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=kmeans --transformer=hmm --num-clusters=8" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=kmeans --transformer=hmm --num-clusters=10" openmoc-assm-1.6.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=birch --transformer=hmm --num-clusters=2" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=birch --transformer=hmm --num-clusters=4" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=birch --transformer=hmm --num-clusters=6" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=birch --transformer=hmm --num-clusters=8" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=birch --transformer=hmm --num-clusters=10" openmoc-assm-1.6.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=dbscan --transformer=hmm --num-clusters=2" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=dbscan --transformer=hmm --num-clusters=4" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=dbscan --transformer=hmm --num-clusters=6" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=dbscan --transformer=hmm --num-clusters=8" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=dbscan --transformer=hmm --num-clusters=10" openmoc-assm-1.6.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=gmm --transformer=hmm --num-clusters=2" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=gmm --transformer=hmm --num-clusters=4" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=gmm --transformer=hmm --num-clusters=6" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=gmm --transformer=hmm --num-clusters=8" openmoc-assm-1.6.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --predictor=gmm --transformer=hmm --num-clusters=10" openmoc-assm-1.6.pbs
