#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export FIRST=$1
export LEN=$2
export NAME=$3

cd ${FIRST}

printf "${FIRST}_dipy_${NAME}.txt\n${LEN}\n${FIRST}_dipy_sorted_plusstrand.wig\n${FIRST}_dipy_sorted_minusstrand.wig\n" | perl ../captureseqtf_cpdplot.pl
