export ILCSOFT=/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01

# -------------------------------------------------------------------- ---

# ---  Use the same compiler and python as used for the installation   ---

# -------------------------------------------------------------------- ---
export PATH=/usr/bin:/usr/bin:${PATH}
export LD_LIBRARY_PATH=/usr/lib64:/usr/lib:/usr/lib:${LD_LIBRARY_PATH}

export CXX=g++
export CC=gcc


#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/lccd/v01-05"


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export COND_DB_DEBUGLOG="/dev/stdout"
export CondDBMySQL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CondDBMySQL/CondDBMySQL_ILC-0-9-6"
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/lcio/v02-12"
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LD_LIBRARY_PATH"


# #--------------------------------------------------------------------------------
# #     ROOT
# #--------------------------------------------------------------------------------
# export ROOTSYS="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/root/6.08.06"
# export PYTHONPATH="$ROOTSYS/lib:$PYTHONPATH"
# export PATH="$ROOTSYS/bin:$PATH"
# export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/CMake/3.6.3/bin:$PATH"


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/ilcutil/v01-05"
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/Marlin/v01-16"
export PATH="$MARLIN/bin:$PATH"
export MARLIN_DLL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinDD4hep/v00-06/lib/libMarlinDD4hep.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DDMarlinPandora/v00-10/lib/libDDMarlinPandora.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinReco/v01-24-01/lib/libMarlinReco.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/PandoraAnalysis/v02-00-00/lib/libPandoraAnalysis.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/LCFIVertex/v00-07-04/lib/libLCFIVertexProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CEDViewer/v01-15/lib/libCEDViewer.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/Overlay/v00-21/lib/libOverlay.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinFastJet/v00-05-01/lib/libMarlinFastJet.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/LCTuple/v01-11/lib/libLCTuple.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinKinfit/v00-06/lib/libMarlinKinfit.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinTrkProcessors/v02-10/lib/libMarlinTrkProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinKinfitProcessors/v00-04/lib/libMarlinKinfitProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/ILDPerformance/v01-06/lib/libILDPerformance.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/Clupatra/v01-03/lib/libClupatra.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/Physsim/v00-04-01/lib/libPhyssim.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/LCFIPlus/v00-06-08/lib/libLCFIPlus.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/FCalClusterer/v01-00/lib/libFCalClusterer.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/ForwardTracking/v01-13/lib/libForwardTracking.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/ConformalTracking/v01-07/lib/libConformalTracking.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/LICH/v00-01/lib/libLICH.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinTPC/v01-04/lib/libMarlinTPC.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/Garlic/v03-01/lib/libGarlic.so:$MARLIN_DLL"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CLHEP/2.3.4.3"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/RAIDA/v01-09"
export PATH="$RAIDA_HOME/bin:$PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/gear/v01-08"
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4hep_ROOT="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DD4hep/v01-07-01"
export DD4hepINSTALL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DD4hep/v01-07-01"
export DD4HEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DD4hep/v01-07-01"
export DD4hep_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DD4hep/v01-07-01"
export PYTHONPATH="$DD4HEP/python:$DD4HEP/DDCore/python:$PYTHONPATH"
export PATH="$DD4HEP/bin:$PATH"
export LD_LIBRARY_PATH="$DD4HEP/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/geant4/10.03.p02"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
export G4SYSTEM="Linux-g++"


#--------------------------------------------------------------------------------
#     QT
#--------------------------------------------------------------------------------
export QTDIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/QT/4.7.4"
export QMAKESPEC="$QTDIR/mkspecs/linux-g++"
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/xercesc/3.1.4"
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/afs/desy.de/project/ilcsoft/sw/boost/1.58.0"


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/KalTest/v02-04"
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     aidaTT
#--------------------------------------------------------------------------------
export AIDATT="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/aidaTT/v00-09"
export PATH="$AIDATT/bin:$PATH"
export LD_LIBRARY_PATH="$AIDATT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/gsl/2.1"
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GBL
#--------------------------------------------------------------------------------
export GBL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/GBL/V02-01-01"
export LD_LIBRARY_PATH="$GBL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinUtil/v01-15/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CED/v01-09-02/bin:$PATH"
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CED/v01-09-02/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/PandoraPFANew/v03-09-00"
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraAnalysis
#--------------------------------------------------------------------------------
export PANDORA_ANALYSIS_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/PandoraAnalysis/v02-00-00"
export PATH="$PANDORA_ANALYSIS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$PANDORA_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/LCFIVertex/v00-07-04/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CEDViewer/v01-15/bin:$PATH"
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/CEDViewer/v01-15/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/FastJet/3.2.1"
export PATH="$FastJet_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$FastJet_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinTPC
#--------------------------------------------------------------------------------
export MARLINTPC="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/MarlinTPC/v01-04"
export PATH="$MARLINTPC/bin:$PATH"


#--------------------------------------------------------------------------------
#     lcgeo
#--------------------------------------------------------------------------------
export lcgeo_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/lcgeo/v00-16-01"
export PYTHONPATH="$lcgeo_DIR/lib/python:$PYTHONPATH"
export PATH="$lcgeo_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$lcgeo_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     DD4hepExamples
#--------------------------------------------------------------------------------
export DD4hepExamples="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/DD4hepExamples/v01-07-01"
export PATH="$DD4hepExamples/bin:$PATH"
export LD_LIBRARY_PATH="$DD4hepExamples/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MySQL
#--------------------------------------------------------------------------------
export MYSQL_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/mysql/usr"
export MYSQL_LIBDIR="$MYSQL_HOME/lib64/mysql"
export MYSQL_PATH="$MYSQL_HOME"
export MYSQL="$MYSQL_HOME"
export PATH="$MYSQL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MYSQL_HOME/lib64:$MYSQL_HOME/lib:$MYSQL_HOME/lib64/mysql:$MYSQL_HOME/lib/mysql:$LD_LIBRARY_PATH"

# --- source GEANT4 INIT script ---
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }
