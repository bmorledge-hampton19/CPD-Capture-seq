#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export ZERO=$1
export NAKEDNA=$2 
export LEN=$3
export SCALE=$4
export NORM=$5
export ORDER=$6
export NAME=$7

cd ${NAKEDNA}
printf "${NAKEDNA}_dipy_${NAME}plot.txt\n${LEN}\n${NAKEDNA}_dipy_sorted_plusstrand.wig\n${NAKEDNA}_dipy_sorted_minusstrand.wig\n" | perl ../capturesequ_clusterplot.pl

cd ..
cd ${ZERO}

printf "${ZERO}_dipy_${NAME}plot.txt\n${LEN}\n${ZERO}_dipy_sorted_plusstrand.wig\n${ZERO}_dipy_sorted_minusstrand.wig\n" | perl ../capturesequ_clusterplot.pl

printf "${ZERO}_dipy_${NAME}plot_${LEN}bp.txt\n../${NAKEDNA}/${NAKEDNA}_dipy_${NAME}plot_${LEN}bp.txt\n${NORM}\n" | perl ../normdiff_humancluster.pl
 
printf "${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff.txt\n${SCALE}\n0\n" | perl ../arbnormcenter_humanclusters.pl

mv ${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff_scal${SCALE}_center0.txt ${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff_scal${SCALE}.txt

printf "../${ORDER}\n${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff_scal${SCALE}.txt\n" | perl ../order_cluster.pl

perl ../formathuman_CDT.pl <${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff_scal${SCALE}_filesort.txt >${ZERO}_dipy_${NAME}plot_${LEN}bp_normdiff_scal${SCALE}_filesort.cdt

