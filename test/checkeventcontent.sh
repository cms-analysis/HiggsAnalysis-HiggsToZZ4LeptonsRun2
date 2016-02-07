#!/bin/sh

source /cmshome/nicola/logincms_amd.sh 
export SCRAM_ARCH=slc5_amd64_gcc462
cd /cmshome/nicola/tmp/test/Moriond/CMSSW_5_3_4/src/HiggsAnalysis/HiggsToZZ4Leptons/test/
eval `scramv1 runtime -sh`

rm -f tmp.txt
rm -f checkeventcontent.txt

cat edmDumpEventContent_a.txt | sed "s?drop ?drop_?g"| sed "s?keep ?keep_?g"  > edmDumpEventContent.txt_tmp

for line in `less edmDumpEventContent.txt_tmp`; do

 echo $line | sed "s?drop_?drop ?g" | sed "s?keep_?keep ?g" >> tmp.txt  

 echo "Line is $line"

 # cat tmp.txt
 cat fragment1_a.txt tmp.txt fragment2.txt  > ../python/hTozzTo4leptons_EventContentReduced_cff.py
 # cmsRun HiggsToZZ_mc_noskim_EScaleCalib_EA2012_Regr_Mor_saveEDM.py
 cmsRun step1_saveEDM.py
 cmsRun HiggsToZZ_mc_noskim_EScaleCalib_EA2012_Regr_Mor.py
 if [ $? != 0 ]; then
   echo "Error: the collection $line cannot be dropped" 
   echo "Error: the collection $line cannot be dropped" >> checkeventcontent.txt
   nlines=`wc -l tmp.txt | awk '{print $1}'`;
   (( m = ${nlines} - 1 ))
   echo "m is $m"
   rm -f tmp2.txt
   cat tmp.txt | tail -n $m > tmp2.txt
   echo $line | sed "s?drop_?keep ?g" >> tmp2.txt
   mv tmp2.txt tmp.txt
   cat fragment1_a.txt tmp.txt fragment2.txt  > ../python/hTozzTo4leptons_EventContentReduced_cff.py
 fi
done
