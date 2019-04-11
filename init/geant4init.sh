#!bash/sh

# . /opt/products/geant4/9.4.p01/env.sh
# 
# #set G4WORKDIR to current directory:
# unset G4WORKDIR
# export G4WORKDIR=$PWD
# export PATH="$G4WORKDIR/bin/Linux-g++:$PATH"
# echo "G4WORKDIR set to $G4WORKDIR"

#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc44_sl6/v01-17-06/geant4/9.5.p02/share/Geant4-9.5.2/geant4make/geant4make.sh

source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc49_sl6/v02-00-01/geant4/10.03.p02/share/Geant4-10.3.2/geant4make/geant4make.sh

export G4WORKDIR=$PWD
echo $PWD
export QTHOME=$QTDIR


 
