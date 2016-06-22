#!/bin/bash


mkdir -p /lustre/cms/store/user/defilip/MonoHiggs/76X/jobdir
mkdir -p /lustre/cms/store/user/defilip/MonoHiggs/76X/histodir

echo "Running HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"
source /cvmfs/cms.cern.ch/cmsset_default.sh

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export PATH=path:$PATH

#exedir=`echo /lustre/home/nicola/slc6/MonoHiggs/Analysis13TeV/Sync13TeV/CMSSW_7_6_3_patch2/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/`
#export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH
#cd ${exedir}
#eval `scramv1 runtime -sh`

export CMSSW_BASE=CMSSW_BASE_DIR
melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc493/
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

#export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
#export PATH=path:$PATH


if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR`;
    cd ${workdir};
else 
    workdir=`echo $PWD`;
    cd ${workdir};
fi

#if [ ! -d "LineShape" ]; then 
# # echo ciao 
# cp -r ${exedir}/LineShape .
#fi

#if [ ! -d "LineShapeNew" ]; then
# # echo ciao
# cp -r ${exedir}/LineShapeNew .
#fi

#if [ ! -d "VBFLineShape" ]; then
# cp -r ${exedir}/VBFLineShape .
#fi

#if [ ! -d "PUvertices" ]; then
# cp -r ${exedir}/PUvertices .
#fi

#if [ ! -d "HqTweights" ]; then
# cp -r ${exedir}/HqTweights .
#fi

#if [ ! -f "CombinedMethod_ScaleFactors_RecoIdIsoSip.root" ]; then 
#	cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip.root .
#fi
#if [ ! -f "CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root" ]; then 
#	cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root .
#fi
#if [ ! -f "MuonScaleFactors_2011_2012_v2.root" ]; then 
#	cp ${exedir}/MuonScaleFactors_2011_2012_v2.root .
#fi
#if [ ! -f "puProfile_Summer12_53X.root"  ]; then 
#	cp ${exedir}/puProfile_Summer12_53X.root .
#fi
#if [ ! -f "pileupMC_2012.root" ]; then 
#	cp ${exedir}/pileupMC_2012.root .
#fi
#if [ ! -f "puProfile_Data_8TeV.root" ]; then 
#	cp ${exedir}/puProfile_Data_8TeV.root .
#fi
#if [ ! -f "ebeOverallCorrections.LegacyPaper.42x.root" ]; then 
#	cp ${exedir}/ebeOverallCorrections.LegacyPaper.42x.root .
#fi
#if [ ! -f "ebeOverallCorrections.Legacy2013.v0.root" ]; then 
#	cp ${exedir}/ebeOverallCorrections.Legacy2013.v0.root .
#fi

# Pileup reweighting
# cp ${exedir}/pileup_MC_Data_76x_50ns_25ns_silver.root .

savedir=`echo /lustre/cms/store/user/defilip/MonoHiggs/76X/histodir`

echo "Working dir is $workdir"
#echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compilereference.sh 4mu


./RunReferenceAnalysis ./sig_input_h150.txt 1 ./bkg_input.txt 1 ./data_input.txt 1 site year mc >& ${workdir}/HZZ4LeptonsAnalysis_log
cp -f ${workdir}/HZZ4LeptonsAnalysis_log /lustre/cms/store/user/defilip/MonoHiggs/76X/jobdir/HZZ4LeptonsAnalysis_log

mv ${workdir}/output.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.
mv ${workdir}/output_txt_vbf.txt ${savedir}/.


