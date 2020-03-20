# PYTHIA configuration file.
# Generated on Fri Mar 20 13:05:24 EDT 2020 with the user supplied options:
# --enable-shared
# --with-boost=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/boost/1.63.0-gnimlf
# --with-hepmc2=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/hepmc/2.06.07-omkpbe2
# --with-lhapdf6=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3
# --with-lhapdf6-plugin=LHAPDF6.h

# Install directory prefixes.
PREFIX_BIN=/data/users/eno/emg_8080/CMSSW_10_2_21/pythia8/bin
PREFIX_INCLUDE=/data/users/eno/emg_8080/CMSSW_10_2_21/pythia8/include
PREFIX_LIB=/data/users/eno/emg_8080/CMSSW_10_2_21/pythia8/lib
PREFIX_SHARE=/data/users/eno/emg_8080/CMSSW_10_2_21/pythia8/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=./
FASTJET3_LIB=./

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/bin/
HEPMC2_INCLUDE=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/include
HEPMC2_LIB=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/hepmc/2.06.07-omkpbe2/lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=true
LHAPDF6_BIN=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/bin/
LHAPDF6_INCLUDE=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/include
LHAPDF6_LIB=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1-gnimlf3/lib

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=
ROOT_INCLUDE=./
ROOT_LIB=./

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=true
BOOST_BIN=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/boost/1.63.0-gnimlf/bin/
BOOST_INCLUDE=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/boost/1.63.0-gnimlf/include
BOOST_LIB=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/boost/1.63.0-gnimlf/lib

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
