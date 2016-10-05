#!/bin/bash

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=agglomerative --num-clusters=2" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=agglomerative --num-clusters=4" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=agglomerative --num-clusters=6" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=agglomerative --num-clusters=8" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=agglomerative --num-clusters=10" openmoc-refl.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=kmeans --num-clusters=2" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=kmeans --num-clusters=4" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=kmeans --num-clusters=6" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=kmeans --num-clusters=8" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=kmeans --num-clusters=10" openmoc-refl.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=birch --num-clusters=2" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=birch --num-clusters=4" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=birch --num-clusters=6" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=birch --num-clusters=8" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=birch --num-clusters=10" openmoc-refl.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=dbscan --num-clusters=2" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=dbscan --num-clusters=4" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=dbscan --num-clusters=6" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=dbscan --num-clusters=8" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=dbscan --num-clusters=10" openmoc-refl.pbs

qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=gmm --num-clusters=2" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=gmm --num-clusters=4" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=gmm --num-clusters=6" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=gmm --num-clusters=8" openmoc-refl.pbs
qsub -P moose -v flags="-a 128 -s 0.05 --num-fine=70 --batchwise=perfect --clusterizer=pinch --transformer=hmm --predictor=gmm --num-clusters=10" openmoc-refl.pbs
