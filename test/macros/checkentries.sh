#!/bin/bash

# for file in `ls  /lustre/cms/store/user/defilip/Summer12_53X_paper_step_analysis_merged/roottree_leptons_QCD*_*8TeV*.root`; do
# for file in `ls  /lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged/roottree_leptons_GluGluToH*7TeV*.root`; do
# for file in `ls  /lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged_local/roottree_leptons_*.root`; do
# for file in `cat filelist_4mu_2012_paper_Bari.txt`; do
# for file in `ls -lstrd1 /lustre/cms/store/user/defilip/histos2e2mu/*.root | grep -v bnn | awk '{print $10}'`; do
for file in `ls  /lustre/cms/store/user/defilip/Summer12_53X_paper_step_analysis_merged_local/roottree_leptons_Summer12_53X_step1_paper_VBFToHToZZTo4L_M-600*`; do
echo $file
cat checkentries.C | sed "s?filename?${file}?g" > tmp.C
g++ -I $ROOTSYS/include tmp.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o checkentries
./checkentries

done

