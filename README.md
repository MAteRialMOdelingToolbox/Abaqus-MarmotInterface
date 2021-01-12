## Abaqus Marmot Interface

Interface files for the usage of Marmot in Abaqus 2019.

# How to use Abaqus together with Marmot

1. Copy one of the environment files (```abaqus_v6.env_gcc_compiler_linux``` is recommended) to your working directory and rename it to ```abaqus_v6.env```.
2. Set environment variable ```MARMOT_DIR``` to your ```Marmot``` directory containing the precompiled Marmot library.
3. Call Abaqus, and provide user.cpp 

# Changes for Abaqus v<2019

To make the interface compatible with older versions of Abaqus, change all ```FOR_NAME``` preprocessor macros to the old specification ( e.g. ```FOR_NAME( umat, UMAT )``` to ```FOR_NAME( umat )``` ).
