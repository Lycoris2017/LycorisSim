#!bash/sh

# -------------------------------------------------------------------- ---
# ---  Use the same compiler and python as used for the installation   ---
# -------------------------------------------------------------------- ---
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin:/cvmfs/sft.cern.ch/lcg/releases/LCG_87/Python/2.7.10/x86_64-slc6-gcc49-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/lib64:/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/lib:/cvmfs/sft.cern.ch/lcg/releases/LCG_87/Python/2.7.10/x86_64-slc6-gcc49-opt/lib:${LD_LIBRARY_PATH}

export CXX=g++
export CC=gcc

#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/lcio/v02-12"
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     ROOT
#--------------------------------------------------------------------------------
export ROOTSYS="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/root/6.08.06"
export PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin:$PATH"

#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/CLHEP/2.3.4.3"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"

#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/gear/v01-08"
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4hep_ROOT="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/DD4hep/v01-07-01"
export DD4hepINSTALL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/DD4hep/v01-07-01"
export DD4HEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/DD4hep/v01-07-01"
export DD4hep_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/DD4hep/v01-07-01"
export PYTHONPATH="$DD4HEP/python:$DD4HEP/DDCore/python:$PYTHONPATH"
export PATH="$DD4HEP/bin:$PATH"
export LD_LIBRARY_PATH="$DD4HEP/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/geant4/10.03.p02"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
export G4SYSTEM="Linux-g++"

source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/geant4/10.03.p02/share/Geant4-10.3.2/geant4make/geant4make.sh

#--------------------------------------------------------------------------------
#     QT
#--------------------------------------------------------------------------------
export QTDIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/QT/4.7.4"
export QMAKESPEC="$QTDIR/mkspecs/linux-g++"
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"

#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/afs/desy.de/project/ilcsoft/sw/boost/1.58.0"


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/xercesc/3.1.4"
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"

