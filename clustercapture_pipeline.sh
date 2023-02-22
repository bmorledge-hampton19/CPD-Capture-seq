#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

export NA=$1
export LEN=$2
export SCALE=$3
export CENTER=$4
cd ${NA}

# set backbground for dipyrimidine reads
printf "${NA}_dipy_orfplot.txt\n${LEN}\n${NA}_dipy_sorted_plusstrand.wig\n${NA}_dipy_sorted_minusstrand.wig\n" | perl ../capturesequ_clusterplot.pl

printf "${NA}_dipy_orfplot_${LEN}bp.txt\n${SCALE}\n${CENTER}\n" | perl ../arbnormcenter_humanclusters.pl

printf "${NA}_dipy_orfplot_${LEN}bp_scal${SCALE}_center${CENTER}.txt" | perl ../cluster_captureseq_order.pl


perl ../formathuman_CDT.pl <${NA}_dipy_orfplot_${LEN}bp_scal${SCALE}_center${CENTER}_order.txt >${NA}_dipy_orfplot_${LEN}bp_scal${SCALE}_center${CENTER}_order.cdt

