#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export ZERO=$1
export NAKEDNA=$2 
export LEN=$3
export NAME=$4

cd ${NAKEDNA}
printf "${NAKEDNA}_dipy_${NAME}.txt\n${LEN}\n${NAKEDNA}_dipy_sorted_plusstrand.wig\n${NAKEDNA}_dipy_sorted_minusstrand.wig\n" | perl ../captureseqtf_cpdplot.pl

cd ..
cd ${ZERO}

printf "${ZERO}_dipy_${NAME}.txt\n${LEN}\n${ZERO}_dipy_sorted_plusstrand.wig\n${ZERO}_dipy_sorted_minusstrand.wig\n" | perl ../captureseqtf_cpdplot.pl

cat ${ZERO}_dipy_${NAME}_${LEN}bp.txt ../${NAKEDNA}/${NAKEDNA}_dipy_${NAME}_${LEN}bp.txt >${ZERO}_dipy_${NAME}_${LEN}bp_normNakeDNA.txt
