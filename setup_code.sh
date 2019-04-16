#!bash/bash


read -p "Enter your OS sys (sl6 or ubuntu): " typeos
echo "Your OS sys is $typeos"


# dimitra's init files:
#source init/setupILCsoft.sh
#source init/geant4init.sh
# mengqing's cvmfs minimum ilcsoft-v02-00-01 init file:
#source init/G4init.sh

if [[ "$typeos" == *"sl"* ]]; then
	echo "Your OS is Science Linux -- $typeos"
	# if you are on a sl6 machine:
	source /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-01/init_ilcsoft.sh
	
elif [[ "$typeos" == *"ubuntu"* ]]; then
	echo "Your OS is Ubuntu -- $typeos"
	# if you are on a ubuntu16.04 machine:
	#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc54_ub1604/v02-00-01/init_ilcsoft.sh
	source init/setupILCsoft_u1604_mod.sh
else
	echo "Input OS is Not valid!"
fi

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-01/init_ilcsoft.sh
