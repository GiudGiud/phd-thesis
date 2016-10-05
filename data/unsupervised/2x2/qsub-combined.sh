#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=agglomerative --num-clusters=2" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=agglomerative --num-clusters=4" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=agglomerative --num-clusters=6" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=agglomerative --num-clusters=8" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=agglomerative --num-clusters=10" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=kmeans --num-clusters=2" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=kmeans --num-clusters=4" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=kmeans --num-clusters=6" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=kmeans --num-clusters=8" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=kmeans --num-clusters=10" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=birch --num-clusters=2" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=birch --num-clusters=4" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=birch --num-clusters=6" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=birch --num-clusters=8" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=birch --num-clusters=10" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=dbscan --num-clusters=2" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=dbscan --num-clusters=4" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=dbscan --num-clusters=6" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=dbscan --num-clusters=8" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=dbscan --num-clusters=10" openmoc-2x2.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=gmm --num-clusters=2" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=gmm --num-clusters=4" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=gmm --num-clusters=6" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=gmm --num-clusters=8" openmoc-2x2.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=combined --transformer=hmm --predictor=gmm --num-clusters=10" openmoc-2x2.pbs
