Interface files for Abaqus 2017 with bft user subroutines (bftUserLibrary)

Instructions:

1. Copy one of the environment files (gcc recommended) to you working directory, and rename to abaqus_v6.env
2. Set environment variable ```BFT_CONSTITUTIVE_MODELLING ``` to your ```constitutiveModelling``` directory containing the precompiled bftUserLibrary directory.
3. Call Abaqus, and provide user.cpp 

To make the interface compatible with older versions of Abaqus, change all ```FOR_NAME``` preprocessor macros to the old specification ( e.g. ```FOR_NAME( umat, UMAT )``` to ```FOR_NAME( umat )``` ).
