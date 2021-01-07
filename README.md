Interface files for Abaqus 2019 with bft user subroutines (Marmot)

Instructions:

1. Copy one of the environment files (gcc recommended!) to you working directory, and rename to abaqus_v6.env
2. Set environment variable ```MARMOT_DIR``` to your ```Marmot``` directory containing the precompiled Marmot library.
3. Call Abaqus, and provide user.cpp 

To make the interface compatible with older versions of Abaqus, change all ```FOR_NAME``` preprocessor macros to the old specification ( e.g. ```FOR_NAME( umat, UMAT )``` to ```FOR_NAME( umat )``` ).
