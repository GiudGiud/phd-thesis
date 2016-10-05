#!/bin/bash

qsub -P moose -v flags="full-core U238 70 10 4" plot-batch-dev-beavrs.pbs
qsub -P moose -v flags="full-core U238 2 10 4" plot-batch-dev-beavrs.pbs
qsub -P moose -v flags="full-core U235 2 10 4" plot-batch-dev-beavrs.pbs