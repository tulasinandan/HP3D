#!/bin/bash
#PBS -N OTV3D
#PBS -A UDEL0004
#PBS -q premium
#PBS -l walltime=6:00:00
#PBS -j oe	
#PBS -S /bin/csh
#PBS -m abe
#PBS -o OTV3D.out
#PBS -e OTV3D.err

ml reset
ml mpt/2.15f

export CODEDIR=$HOME/HP3D
export PAR=OTV3D  # paramfile = param_$PAR
export NSTEPS=10
export RESTART=0 # Comment this out to start from scratch.
export FILENO=0


export PROGNAME='p3d'

########## END OF USER MODIFIED VARIABLES  ####################
###############################################################

##########  Various initial things #####################

export DR=$SCRATCH/$PAR
export NPES=`\echo 'pex*pey*pez/36 + 1' | cat param_$PAR - | cpp | tail -1 | bc`
export NPROCS=`\echo 'pex*pey*pez' | cat param_$PAR - | cpp | tail -1 | bc`

############ Initialize RUN ####################
 
if [ "$1" == "init" ] ; then
 case $RESTART in 
  0) #### If RESTART=0
   if [ -e $DR ]; then echo 'Directory Exists'; exit; fi
   mkdir -p $DR
 
 #  ########## Copy files into scratch directory and compile
#set -x
 
    cp $0 $DR
    cd $DR ; mkdir compiled; mkdir FPData BData
    cp $CODEDIR/runtime/param_$PAR .;
    cp $CODEDIR/src/*.h $CODEDIR/src/*.M \
       $CODEDIR/src/*.F90 $CODEDIR/src/makefile \
       compiled/.
  
    #*/###########################################
    cd $DR/compiled # Change to compile directory
 
    # cp param file to simply "param" so cpp will include it. 
    cp ../param_$PAR param
 
    # Untar program files
 #  tar -xvf ${VERSION}.tar
 
    # Figure out what kind of machine I'm running on. Right now, this 
    # script is only set up to run on the IBM. 
    if [ "`uname -s`" = "AIX" ]; then export MACH='AIX'
    elif [ "`uname -s`" = "Linux" ]; then export MACH='Linux'
    fi
 
    # Determine which initialization routine is being used. 
    export INITFILE=`awk '/^#define/ && / init_scheme / {print $3}' param`
 
   echo " "
   echo Compiling the source
   echo " "
   echo " "

   # Make init program which outputs initial dump file. 
   make par=$PAR mach=$MACH initfile=$INITFILE init-hybrid
   # Make p3d
   make par=$PAR mach=$MACH initfile=$INITFILE p3d-hybrid
   # Make moviecombine
   make par=$PAR mach=$MACH initfile=$INITFILE moviecombine
   # Make moviecombinefull
   make par=$PAR mach=$MACH initfile=$INITFILE moviecombinefull
   # Make vidst2d
   make par=$PAR mach=$MACH vdist2d

    mv consistency_check init-hybrid p3d-hybrid moviecombine moviecombinefull vdist2d ../ #move executables from compile dir.
    \rm *.o *.mod
  esac
 
cd $DR
qsub -l select=${NPES}:ncpus=36:mpiprocs=36 $0
fi
 ##########  End of initialization ####################

#############################################################
###########  Begin Run  #####################################
 
 if [ "$1" == "" ] ; then
 
 set -x
 echo $0
    cd $DR
 
    echo "#################### RUNNING  PROGRAM ####################"
    if [ "$RESTART" = "0" ]; then 
    mpiexec_mpt -n $NPROCS omplace $DR/init-hybrid >& init.stdout
     echo "####################### CREATED THE INITIAL FILES ######################"
    fi
    mpiexec_mpt -n $NPROCS omplace $DR/p3d-hybrid $FILENO $NSTEPS >& p3d.stdout
    echo "#################### FINISHED PROGRAM ####################"
 
    ln -s $DR/param_$PAR $DR/paramfile
 
    # Rename movie files so that they are compatible with Shay's idl codes. 
    /bin/cp $DR/param_${PAR} $DR/BData/paramfile $DR/FPData/paramfile
#   cd $DR/BData; time ../moviecombine $NUMSLICES 1 $NUMSLICES 1 >>comb_out
#   for movvar in bx by bz jix jiy jiz log jex jey jez pe pxx pyy pzz pxy pxz pyz \
#           jtotx jtoty jtotz n
#     do
#     if [ -e $movvar ]; then
#        mv -f $movvar movie.${movvar}_${RN}.in
#     fi
#     done

#   cd $DR/FPData; time ../moviecombinefull $NUMSLICESFULL 1 $NUMSLICESFULL 1 >> comb_out

    cd $DR; mkdir dmp; mv p3d-*.* dmp
    
    exit 0
 
 fi
