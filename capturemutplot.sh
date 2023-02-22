#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export FIRST=$1
export NAME=$2

cd ${FIRST}

printf "${FIRST}_dipy_${NAME}.txt\n${FIRST}_dipy_sorted_plusstrand.wig\n${FIRST}_dipy_sorted_minusstrand.wig\n" | perl ../captureseq_mutplot_indivpos.pl

perl ../expand_recurrent_ets.pl <${FIRST}_dipy_${NAME}.bed >${FIRST}_dipy_${NAME}_expand.bed

fastaFromBed -s -nameOnly -fi ../hg19_puc19.fa -bed ${FIRST}_dipy_${NAME}_expand.bed -fo ${FIRST}_dipy_${NAME}_expand.fa

perl ../analyze_recurrent_ets.pl <${FIRST}_dipy_${NAME}_expand.fa >${FIRST}_dipy_${NAME}_expand_ets.txt


