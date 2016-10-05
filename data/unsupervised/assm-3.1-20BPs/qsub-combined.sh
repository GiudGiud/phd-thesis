#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=agglomerative --transformer=hmm --num-clusters=2" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=agglomerative --transformer=hmm --num-clusters=4" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=agglomerative --transformer=hmm --num-clusters=6" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=agglomerative --transformer=hmm --num-clusters=8" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=agglomerative --transformer=hmm --num-clusters=10" openmoc-assm-3.1-20BPs.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=kmeans --transformer=hmm --num-clusters=2" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=kmeans --transformer=hmm --num-clusters=4" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=kmeans --transformer=hmm --num-clusters=6" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=kmeans --transformer=hmm --num-clusters=8" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=kmeans --transformer=hmm --num-clusters=10" openmoc-assm-3.1-20BPs.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=birch --transformer=hmm --num-clusters=2" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=birch --transformer=hmm --num-clusters=4" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=birch --transformer=hmm --num-clusters=6" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=birch --transformer=hmm --num-clusters=8" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=birch --transformer=hmm --num-clusters=10" openmoc-assm-3.1-20BPs.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=dbscan --transformer=hmm --num-clusters=2" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=dbscan --transformer=hmm --num-clusters=4" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=dbscan --transformer=hmm --num-clusters=6" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=dbscan --transformer=hmm --num-clusters=8" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=dbscan --transformer=hmm --num-clusters=10" openmoc-assm-3.1-20BPs.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=gmm --transformer=hmm --num-clusters=2" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=gmm --transformer=hmm --num-clusters=4" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=gmm --transformer=hmm --num-clusters=6" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=gmm --transformer=hmm --num-clusters=8" openmoc-assm-3.1-20BPs.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --predictor=gmm --transformer=hmm --num-clusters=10" openmoc-assm-3.1-20BPs.pbs
