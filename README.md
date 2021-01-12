# Abaqus Marmot Interface

This interface files enable the usage of [Marmot](https://github.com/MAteRialMOdelingToolbox/marmot) in the FEA software package Abaqus.


## How to use Abaqus together with Marmot

1. Copy one of the environment files (```abaqus_v6.env_gcc_compiler_linux``` is recommended) to your working directory and rename it to ```abaqus_v6.env```.
2. Set environment variable ```MARMOT_DIR``` to your ```Marmot``` directory containing the precompiled ```Marmot``` library. The default value of ```MARMOT_DIR``` will be set to ```/usr/local```.
3. Call your Abaqus job file via the command line ```abaqus j=myJob.inp user=user.cpp interactive```. 

## Necessary changes for the usage with Abaqus version<2019

To make the interface compatible with older versions of Abaqus, change all ```FOR_NAME``` preprocessor macros to the old specification ( e.g. ```FOR_NAME( umat, UMAT )``` to ```FOR_NAME( umat )``` ).

## How to define material properties for using a material model of Marmot

A Marmot material is defined as follows in your ```.inp```-file. For the sake of simplicity, a simple ```LinearElastic``` material model from Marmot is used.

```abaqus
**                                       _   
**  _ __ ___   __ _ _ __ _ __ ___   ___ | |_ 
** | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
** | | | | | | (_| | |  | | | | | | (_) | |_ 
** |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
**
** 
** Define your desired user material; The following lines may be simply copied 
** and modified for your purposes.
** 
** --------------------------------------------------------------------------
** 
** 
** name ... name of your user-material (not case sensitive, since the Abaqus 
**          input reader is not); if you have multiple set of properties just
**          add -1, -2 as a suffix, e.g., LinearElastic-1, LinearElastic-2
** 
*Material, name=LinearElastic-1
** 
** Define the number of state variables; The number of required state variables 
** can be found in the header file of your material model defined within the 
** method getNumberOfRequiredStateVars(). If a too small number is chosen, an 
** exception is thrown anyhow.
** 
*Depvar
1
** 
** The following lines are optional and enable the description of the individual
** state variables.
** 
1, eps11, "normal strain in horizontal direction"
** 
** In the following lines, the material properties of the user material are 
** defined. If the system matrix will be unsymmetric (e.g. for non-associated
** plasticity models), the unsymm keyword must be added. To determine the number 
** of material properties the parameter constants must be specified. The required 
** number of material properties can be found for each Marmot material in the 
** constructor. An assert is called in Marmot, when the number of constants is
** chosen too small.
**
*User Material, unsymm, Constants=2
** 
** In the following lines, the material properties are specified. Please make 
** sure that not more than 8 (!) numbers are specified per line (this is a 
** convention from Abaqus). If you need more than 8 parameters, start with a
** new line.
**
** E    nu
30000, 0.2
```

## How to define element properties for using a finite element defined in Marmot

todo


