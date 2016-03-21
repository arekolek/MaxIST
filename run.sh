#!/bin/bash

# TODO save random seed

if [ -z $OMP_NUM_THREADS ]
then
  OMP_NUM_THREADS=12
fi

echo -n 'Question prompting the experiment: '
RATIONALE=$(head -n 1 | awk '{print "# " $0}')

JOB_ID=$(env TS_ONFINISH=./run_after.sh ts -E -N $OMP_NUM_THREADS ./eval.exe $@)

echo "$RATIONALE
# 
# Commit: $(git rev-parse HEAD)" > output/ts_$JOB_ID.txt

