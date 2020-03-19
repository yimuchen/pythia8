#!/bin/bash

installPythia() {
	CORES=$1
	if [ -z ${CORES} ]; then
		CORES=1
	fi

	if [ -z "$CMSSW_BASE" ]; then
		echo "$CMSSW_BASE required"
		exit 1
	fi

	cd ${CMSSW_BASE}
	# some really bad ways to get info out of scram
	HEPMC_BASE=$(scram tool info hepmc | grep "HEPMC_BASE" | sed 's/HEPMC_BASE=//')
	BOOST_BASE=$(scram tool info boost | grep "BOOST_BASE" | sed 's/BOOST_BASE=//')
	LHAPDF_BASE=$(scram tool info lhapdf | grep "LHAPDF_BASE" | sed 's/LHAPDF_BASE=//')

	# get pythia8 source and compile
	git clone git@github.com:kpedro88/pythia8 -b emg/230
	cd pythia8
	# configure for c++11 if 7_1_X
	EXTRA=""
	case $CMSSW_VERSION in
	CMSSW_7_1_*)
		EXTRA='--cxx-common="-std=c++11 -fPIC"'
	;;
	esac
	./configure --enable-shared --with-boost=${BOOST_BASE} --with-hepmc2=${HEPMC_BASE} --with-lhapdf6=${LHAPDF_BASE} --with-lhapdf6-plugin=LHAPDF6.h "$EXTRA"
	make -j ${CORES}
	make install

	# create xml for tool
	cd $CMSSW_BASE
	cat << 'EOF_TOOLFILE' > pythia8.xml
<tool name="pythia8" version="230-emg">
  <lib name="pythia8"/>
  <client>
    <environment name="PYTHIA8_BASE" default="$CMSSW_BASE/pythia8"/>
    <environment name="LIBDIR" default="$PYTHIA8_BASE/lib"/>
    <environment name="INCLUDE" default="$PYTHIA8_BASE/include"/>
  </client>
  <runtime name="PYTHIA8DATA" value="$PYTHIA8_BASE/share/Pythia8/xmldoc"/>
  <use name="hepmc"/>
  <use name="lhapdf"/>
</tool>
EOF_TOOLFILE

	# install tool in scram
	mv ${CMSSW_BASE}/pythia8.xml ${CMSSW_BASE}/config/toolbox/${SCRAM_ARCH}/tools/selected
	scram setup pythia8

	# update CMSSW dependencies
	case $CMSSW_VERSION in
	CMSSW_7_1_*)
		# better to use 'scram b checkdeps' but this is not available in 71X
		scram b echo_pythia8_USED_BY | tr ' ' '\n' | grep "self" | cut -d'/' -f2-3 | sort -u > pkgs.txt
		cd $CMSSW_BASE/src
		git cms-addpkg -f ../pkgs.txt
	;;
	*)
		cd $CMSSW_BASE/src
		scram b checkdeps
	;;
	esac

	scram b -j ${CORES}
}

installPythia $1
