#!/bin/bash
PREFIX=$1;
ASM_DIR=$2;
READLEN=$3;
READ_SR=$4;
gatekeeper -dumpfragments -tabular ${ASM_DIR}/${PREFIX}.gkpStore |awk '{print $1}' > ${PREFIX}.uid
tigStore -g ${ASM_DIR}/${PREFIX}.gkpStore -t ${ASM_DIR}/${PREFIX}.tigStore 2 -U -d layout > unitig_layout.txt
makeAdjustmentFactorsForNumReadsForAStatBasedOnGC unitig_layout.txt > unitigCoverageAdjustedForGC.txt
cat unitig_layout.txt | compute_sr_cov.revisedForGCContig.pl ${PREFIX}.uid $READ_SR $READLEN superReadSequences_shr.frg unitigCoverageAdjustedForGC.txt 1> unitig_cov.txt 2>global_arrival_rate.txt
tigStore -g ${ASM_DIR}/${PREFIX}.gkpStore -t ${ASM_DIR}/${PREFIX}.tigStore 2 -E unitig_cov.txt 1> tigStore.err 2>&1
tigStore -g ${ASM_DIR}/${PREFIX}.gkpStore -t ${ASM_DIR}/${PREFIX}.tigStore 1 -E unitig_cov.txt 1> tigStore.err 2>&1

