#!/bin/bash


mkdir -p $_CONDOR_SCRATCH_DIR/histodir
mkdir -p /eos/uscms/store/user/ndefilip/jobdir
mkdir -p /eos/uscms/store/user/ndefilip/histodir

workdir=${PWD}
echo "Running HtoZZto4Leptons Analysis with executables RunHZZ4LeptonsAnalysis"
source /uscmst1/prod/sw/cms/bashrc prod
exedir=`echo /uscms/home/ndefilip/nobackup/CMSSW_5_3_9/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/FSR2012/INCLUSIVE_MACROS`
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH
cd ${exedir}
eval `scramv1 runtime -sh`

melalibdir=${CMSSW_BASE}/lib/slc5_amd64_gcc462/
export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH

if [ -d "$_CONDOR_SCRATCH_DIR/" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR/`;
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

if [ ! -f "CombinedMethod_ScaleFactors_RecoIdIsoSip.root" ]; then 
	cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip.root .
fi
if [ ! -f "CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root" ]; then 
	cp ${exedir}/CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root .
fi
if [ ! -f "MuonScaleFactors_2011_2012_v2.root" ]; then 
	cp ${exedir}/MuonScaleFactors_2011_2012_v2.root .
fi
if [ ! -f "puProfile_Summer12_53X.root"  ]; then 
	cp ${exedir}/puProfile_Summer12_53X.root .
fi
if [ ! -f "pileupMC_2012.root" ]; then 
	cp ${exedir}/pileupMC_2012.root .
fi
if [ ! -f "puProfile_Data_8TeV.root" ]; then 
	cp ${exedir}/puProfile_Data_8TeV.root .
fi
if [ ! -f "ebeOverallCorrections.LegacyPaper.42x.root" ]; then 
	cp ${exedir}/ebeOverallCorrections.LegacyPaper.42x.root .
fi
if [ ! -f "ebeOverallCorrections.Legacy2013.v0.root" ]; then 
	cp ${exedir}/ebeOverallCorrections.Legacy2013.v0.root .
fi

# savedir=`echo ${exedir}/log`
savedir=`echo /eos/uscms/store/user/ndefilip/histodir`

echo "Working dir is $workdir"
echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

${exedir}/RunReferenceAnalysis ${exedir}/sig_input_h150.txt 1 ${exedir}/bkg_input.txt 1 ${exedir}/data_input.txt 1 site year mc >& ${workdir}/HZZ4LeptonsAnalysis_log
mv ${workdir}/HZZ4LeptonsAnalysis_log /eos/uscms/store/user/ndefilip/jobdir/HZZ4LeptonsAnalysis_log
# rm -f  ${workdir}/HZZ4LeptonsAnalysis_log /eos/uscms/store/user/ndefilip/jobdir/HZZ4LeptonsAnalysis_log

mv ${workdir}/output.root    ${savedir}/.
mv ${workdir}/output_bnn.txt ${savedir}/.
mv ${workdir}/output_bnn.root ${savedir}/.
mv ${workdir}/output_txt.txt ${savedir}/.
mv ${workdir}/output_txt_vbf.txt ${savedir}/.


