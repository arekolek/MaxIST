#!/bin/bash

FILE=output/ts_$1.txt

ts -i $1 | awk '{print "# " $0}' >> $FILE
cat $3 >> $FILE

