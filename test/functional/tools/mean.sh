#!/bin/bash

DATA_STD=`cat stat.txt | grep ${1} | grep Standard | awk '{print $9}'`
DATA_VAR=`cat stat.txt | grep ${1} | grep Variance | awk '{print $8}'`

NUM=`echo "$DATA_STD" | wc -l`
SUM=`echo $DATA_STD | tr ' ' '+' | bc -l`

NUM_VAR=`echo "$DATA_VAR" | wc -l`
SUM_VAR=`echo $DATA_VAR | tr ' ' '+' | bc -l`

MEAN=`echo "$SUM / $NUM" | bc -l`
VARIANCE=`echo "$SUM_VAR / $NUM_VAR" | bc -l`

echo $MEAN $VARIANCE
