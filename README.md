Interface files for Abaqus with bft user subroutines (bftUserLibrary)

Instructions:

1. Copy one of the environment files (gcc recommended) to you working directory, and rename to abaqus_v6.env
2. Set environment variable ```BFT_CONSTITUTIVE_MODELLING ``` to your ```constitutiveModelling``` directory containing the precompiled bftUserLibrary directory.
3. Call Abaqus, and provide either user.cpp (Abaqus 6.14) or user2017.cpp (Abaqus >= 2016)

