#!/bin/bash

echo "Processing on " `hostname` "at " `date` 

mkdir -p /home/tmp/defilip/$$
mkdir -p /lustre/cms/store/user/defilip/MonoHiggs/76X/jobdir
mkdir -p /lustre/cms/store/user/defilip/MonoHiggs/76X/histodir

workdir=${PWD}
echo "Running HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"
source /cmshome/nicola/slc6/logincms_cvmfs_slc6.sh
export SCRAM_ARCH=slc6_amd64_gcc491
exedir=`echo /cmshome/nicola/slc6/MonoHiggs/Analysis13TeV/Sync13TeV/CMSSW_7_6_3_patch2/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros`
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH
cd ${exedir}
eval `scramv1 runtime -sh`

melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc493/
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

if [ -d "/home/tmp/defilip/$$" ]; then
    workdir=`echo /home/tmp/defilip/$$`;
    cd ${workdir};
fi

if [ ! -d "LineShape" ]; then 
 # echo ciao 
 cp -r ${exedir}/LineShape .
fi

if [ ! -d "LineShapeNew" ]; then
 # echo ciao
 cp -r ${exedir}/LineShapeNew .
fi

if [ ! -d "VBFLineShape" ]; then
 cp -r ${exedir}/VBFLineShape .
fi

if [ ! -d "PUvertices" ]; then
 cp -r ${exedir}/PUvertices .
fi

if [ ! -d "HqTweights" ]; then
 cp -r ${exedir}/HqTweights .
fi

# cp ${exedir}/check_DoubleElectron_2011.overlap .
# cp ${exedir}/check_DoubleElectronMuEG_2012.overlap .
# cp ${exedir}/Electron_scale_factors_IDISOSIP_combined.root .
cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip.root .
cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root .
# cp ${exedir}/SF2013_V5_old2011.root .
# cp ${exedir}/muonid_hcp-05.10.2012-with-tk-v2.root . 
cp ${exedir}/MuonScaleFactors_2011_2012_v2.root .
cp ${exedir}/puProfile_Summer12_53X.root .
cp ${exedir}/pileupMC_2012.root .
cp ${exedir}/puProfile_Data_8TeV.root .
cp ${exedir}/ebeOverallCorrections.LegacyPaper.42x.root .
cp ${exedir}/ebeOverallCorrections.Legacy2013.v0.root .

cp ${exedir}/pileup_MC_Data_76x_50ns_25ns_silver.root .

savedir=`echo /lustre/cms/store/user/defilip/MonoHiggs/76X/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunReferenceAnalysis ${exedir}/sig_input_h150.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 Bari year mc >& ${workdir}/HZZ4LeptonsAnalysis_log
cp -f ${workdir}/HZZ4LeptonsAnalysis_log /lustre/cms/store/user/defilip/MonoHiggs/76X/jobdir/HZZ4LeptonsAnalysis_log
cp -f ${workdir}/output.root    ${savedir}/.
cp -f ${workdir}/output_bnn.txt ${savedir}/.
cp -f ${workdir}/output_bnn.root ${savedir}/.
cp -f ${workdir}/output_txt.txt ${savedir}/.
cp -f ${workdir}/output_txt_vbf.txt ${savedir}/.

# cleaning the worker node
if [ -d "/home/tmp/defilip/$$" ]; then
    rm -f -R *
    rm -f *
fi
