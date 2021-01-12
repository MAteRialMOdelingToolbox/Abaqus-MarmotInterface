# Abaqus Marmot Interface

Interface files for the usage of [Marmot](https://github.com/MAteRialMOdelingToolbox/marmot) in Abaqus 2019.

## How to use Abaqus together with Marmot

1. Copy one of the environment files (```abaqus_v6.env_gcc_compiler_linux``` is recommended) to your working directory and rename it to ```abaqus_v6.env```.
2. Set environment variable ```MARMOT_DIR``` to your ```Marmot``` directory containing the precompiled Marmot library.
3. Call your Abaqus with job via the command line ```abaqus j=myJob.inp cpus=1 user=user.cpp interactive```. 

## How to define material properties for using a material model of Marmot

todo

## How to define element properties for using a finite element defined in Marmot

todo

## Changes for the usage with Abaqus version<2019

To make the interface compatible with older versions of Abaqus, change all ```FOR_NAME``` preprocessor macros to the old specification ( e.g. ```FOR_NAME( umat, UMAT )``` to ```FOR_NAME( umat )``` ).
