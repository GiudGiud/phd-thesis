#!/bin/bash

qsub -P moose -v flags="assm-1.6 U238 70 0 1" plot-batch-var.pbs
qsub -P moose -v flags="assm-1.6 U238 2 0 1" plot-batch-var.pbs
qsub -P moose -v flags="assm-1.6 U235 2 0 1" plot-batch-var.pbs