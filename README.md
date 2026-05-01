# Abaqus Marmot Interface

These interface files enable the usage of [Marmot](https://github.com/MAteRialMOdelingToolbox/marmot) in the FEA software package Abaqus.

## How to use Abaqus together with Marmot

1. Copy one of the environment files (e.g., `abaqus_v6.env_gcc_compiler_linux` is recommended) to your working directory and rename it to `abaqus_v6.env`.
2. Set the environment variable `MARMOT_DIR` to your `Marmot` directory containing the precompiled `Marmot` library. The default value of `MARMOT_DIR` will be set to `/usr/local`.
3. Call your Abaqus job file via the command line: `abaqus j=myJob.inp user=user.cpp interactive`. 

## Necessary changes for usage with Abaqus version < 2019

To make the interface compatible with older versions of Abaqus, change all `FOR_NAME` preprocessor macros to the old specification (e.g. change `FOR_NAME( umat, UMAT )` to `FOR_NAME( umat )`).

---

## How to define a User Material (UMAT) from Marmot

A Marmot material (when used with standard Abaqus elements) is defined as follows in your `.inp`-file. For the sake of simplicity, a `LinearElastic` material model from Marmot is used here.
```abaqus
** name ... name of your user-material (not case sensitive, since the Abaqus 
**          input reader is not); if you have multiple sets of properties just
**          add -1, -2 as a suffix, e.g., LinearElastic-1, LinearElastic-2
** 
*Material, name=LinearElastic-1
** 
** Define the number of state variables; The number of required state variables 
** can be found in the header file of your material model defined within the 
** method getNumberOfRequiredStateVars(). An exception is thrown if this is too low.
** 
*Depvar
1
** 
** (Optional) Describe the individual state variables.
** 
1, eps11, "normal strain in horizontal direction"
** 
** Define the number of material properties. If the system matrix is unsymmetric 
** (e.g. for non-associated plasticity models), the unsymm keyword must be added.
**
*User Material, unsymm, Constants=2
** 
** Specify the material properties. 
** WARNING: Abaqus strictly requires a maximum of 8 numbers per line!
** If you need more than 8 parameters, continue on the next line.
**
** E      nu
25850, 0.18
```

---

## How to define a User Element (UEL) from Marmot

Using a Marmot finite element requires a specific workflow because the C++ interface needs to map Abaqus integer codes to Marmot string names. This is achieved using Abaqus Parameter Tables.

### 1. Define the Parameter Tables
At the model level, define a `TABLE COLLECTION` that maps your arbitrary integer codes to Marmot's internal C++ class names and property counts. 
```abaqus
** Define the schema for our mapping table: ID (int), Name (string), nProps (int)
*PARAMETER TABLE TYPE, NAME=int_string_int, PARAMETERS=3
INTEGER                    
STRING                                                                       
INTEGER

** Create the collection that the C++ interface will read
*TABLE COLLECTION, NAME=UEL_CODES

** Map Element Codes -> Marmot Element Names
*PARAMETER TABLE, TYPE=int_string_int, LABEL=UEL_ELEMENTS
** el-code, Marmot el-name, number of element properties
100,        C3D20R,         0

** Map Material Codes -> Marmot Material Names
*PARAMETER TABLE, TYPE=int_string_int, LABEL=UEL_MATERIALS
** mat-code, Marmot mat-name, number of material properties
999,         VonMises,        6 
```

### 2. Define the User Element (`*USER ELEMENT`)
This block defines the topology and variable sizes for the UEL. Note that it does *not* assign actual material values yet.

```abaqus
** TYPE: The name must be 'U' followed by a number (e.g., U020).
** VARIABLES: Must equal the total number of state variables across all integration points.
**            Example: 8 IPs * (1 GCDP sdv + 6 stress + 6 strain) = 104 in total.
** I PROPERTIES: Must be at least 3 (Element Code, Material Code, Init Flags).
*USER ELEMENT, UNSYMM, TYPE=U020, NODES=20, COORDINATES=3, PROPERTIES=8, VARIABLES=104, I PROPERTIES=3
1,2,3
2,1,2,3
3,1,2,3
4,1,2,3
5,1,2,3
6,1,2,3
7,1,2,3
8,1,2,3
9,1,2,3
10,1,2,3
11,1,2,3
12,1,2,3
13,1,2,3
14,1,2,3
15,1,2,3
16,1,2,3
17,1,2,3
18,1,2,3
19,1,2,3
20,1,2,3
```

### 3. Assign the UEL Properties (`*UEL PROPERTY`)
Finally, assign the properties to an element set. Real properties are defined first, followed by the integer properties.

**Rule for Real Properties:** The material properties are listed first, immediately followed by the element properties (if any). You must adhere to the Abaqus limit of exactly **8 parameters per line**. Pad with zeros or `UNUSED` if necessary to format your lines cleanly.

**Rule for Integer Properties:** The three required integer properties are:
1. Element Code (matches `UEL_ELEMENTS` table)
2. Material Code (matches `UEL_MATERIALS` table)
3. Additional Definitions Flag (Bitwise flag: `0` = None, `1` = Geostatic Stress, `2` = Marmot Initialization, `3` = Both).
```abaqus
*UEL PROPERTY, ELSET=concrete                                                
** ALWAYS EXACTLY 8 PARAMETERS PER LINE!
** First are the material properties (e.g., GCDP), then the element properties.
**  1,      2,      3,      4,      5,          6,      7,          8,  
**  E       nu      fcy     Hlin,   dExpHard,   delta,  UNUSED,     UNUSED,
25850,      0.18,   10.33,  1000,   0,          0,      0,          0
**
** INTEGER PROPERTIES come always after the real (float) properties
** el-code, mat-code, additional definitions flag
100,        999,      0
```
