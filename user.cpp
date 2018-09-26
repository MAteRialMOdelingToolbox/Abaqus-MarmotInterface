#include <aba_for_c.h>
#include "userLibrary.h"
#include <iostream>
#include <string>
#include <sstream>
#include "bftUel.h"
#include "bftMaterialHypoElastic.h"

//These functions are provided to the 'sub' UMATs for easy printing Messages and Warnings. 
//

namespace MainConstants
{
	bool printWarnings = false;    
	bool printMessages = false;

    enum AdditionalDefinitions {
        GeostaticStressDefiniton = 0x01 << 0,
    };
}

class MakeString
{
    public:
        std::stringstream stream;
        operator std::string() const { return stream.str(); }

        template<class T>
        MakeString& operator<<(T const& VAR) { stream << VAR; return *this; }
};

extern "C" 
{
        void stdb_abqerr_(const int *lop, const char* stringZT, const int *intArray, const double *realArray, const char *appendix, const int lengthString, const int lengthAppendix);
        void xit_();
}

extern "C" bool warningToMSG(const std::string& message)
{
    // return always false(!)
    const int lop = -1;
    if(MainConstants::printWarnings)
            stdb_abqerr_(&lop, message.c_str(), nullptr , nullptr , nullptr, message.length(), 0);
    return false;
}

extern "C" bool notificationToMSG(const std::string& message)
{
    // return always true(!)
    const int lop = 1;
    if(MainConstants::printMessages)
			stdb_abqerr_(&lop, message.c_str(), nullptr , nullptr , nullptr, message.length(), 0);
    return true;
}


extern "C" void uel_(
        double rightHandSide[/*lVarx , nRightHandSide*/],           // right hand side load vector(s) 1: common, 2: additional for RIKS (see documentation)
        double KMatrix[/*nDof * nDof*/],                            // stiffness matrix 
        double stateVars[],                                         // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        double energies[8],                                         // may be updated: energies
        const int &nDegreesOfFreedom ,
        const int &nRightHandSide,                
        const int &nStateVars, 
        const double properties[/*nProperties*/],
        const int &nProperties,
        const double coordinates[/*mcrd, nNodes*/],                    // undeformed coordinates of the node respective DOFs
        const int &maxNCoords,                                        // max number of coordinates (see documentation)
        const int &nNodes,                  
        const double U[/*nDof*/],                                   // current solution (end of increment)
        const double dU[/*mlvarx=nDof(?), nRightHandSide*/],        // increment of solutions
        const double UDot[/*nDof*/],                                // first derivative (velocity..)
        const double UDotDot[/*nDof*/],                             // second derivative (acceleration..)
        const int &elementType,                                     // user defined element type id
        const double time[2],                                       // 1: time of step, 2: total time
        const double &dTime,                                        // time increment
        const int& stepNumber,
        const int& incrementNumber,
        const int& elementNumber,                  
        const double solutionParams[],                              // solution procedure dependent parameters 
        const int& nDLoadActive,                                    // id of load / flux currently active on this element
        const int distributedLoadTypes[/*mDLoads, * */],            // An array containing the integers used to define distributed load types for the element.  
        const double distributedLoadMags[/*mDloads, * */],          // magnitudes @ end of increment
        const double Predef[/*nPredef*/], 
        const int &nPredef,
        const int lFlags[],
        const int mlvarx[],
        const double dDistributedLoadMags[/*mDloads, * */],         // increment of magnitudes 
        const int &mDload,                                          // total number of distributed loads and fluxes defined on this element
        double &pNewDT,         
        const int integerProperties[],
        const int &nIntegerProperties, 
        const double &period)
{    
        if ( nIntegerProperties != 5 )
            throw std::invalid_argument( MakeString() << "Uel: insufficient integer properties (" << nIntegerProperties <<") provided, but 5 are required");

        userLibrary::ElementCode elementCode =  static_cast<userLibrary::ElementCode> ( integerProperties[0]);
        userLibrary::MaterialCode materialID =  static_cast<userLibrary::MaterialCode>( integerProperties[1] );
        const int nPropertiesElement =          integerProperties[2];
        const int nPropertiesUmat =             integerProperties[3];
        const int additionalDefinitions =       integerProperties[4];

        const double* propertiesUmat =    &properties[0];
        const double* propertiesElement = &properties[nPropertiesUmat];

        BftUel* myUel = userLibrary::UelFactory(elementCode, propertiesElement, nPropertiesElement, elementNumber, materialID, propertiesUmat, nPropertiesUmat); 

        const int nNecessaryStateVars = myUel->getNumberOfRequiredStateVars();

        if ( nNecessaryStateVars > nStateVars )
            throw std::invalid_argument( MakeString() << "Uel with code " << elementCode << " and material " << materialID << ": insufficient stateVars (" << nStateVars <<") provided, but "<< nNecessaryStateVars <<" are required");

        myUel->assignStateVars(stateVars, nStateVars);

        myUel->initializeYourself(coordinates);

        int additionalDefinitionProperties = 0;
        if( additionalDefinitions & MainConstants::AdditionalDefinitions::GeostaticStressDefiniton ) {
            if(lFlags[0] == Abaqus::UelFlags1::GeostaticStress){
                const double* geostaticProperties = &propertiesElement[ nPropertiesElement + additionalDefinitionProperties ];
                myUel->setInitialConditions( BftUel::GeostaticStress, geostaticProperties ); }
            additionalDefinitionProperties += 5;  
        }

        // compute K and P 
        myUel->computeYourself(U , dU, rightHandSide, KMatrix, time, dTime, pNewDT); 

        // compute distributed loads in nodal forces and add it to P 
        //for (int i =0; i<mDload; i++){
            //if (distributedLoadMags[i]<1.e-16)
                //continue;
            //myUel->computeDistributedLoad(BftUel::Pressure, rightHandSide, distributedLoadTypes[i], &distributedLoadMags[i], time, dTime);}

        delete myUel;
}

extern "C" void umat_(
        /*to be def.*/  double stress[],                // stress vector in order: S11, S22, (S33), S12, (S13), (S23) 
        /*to be def.*/  double stateVars[],        // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        /*to be def.*/  double dStressDDStrain[],  // material Jacobian matrix ddSigma/ddEpsilon
        /*to be def.*/  double &sSE,                    // specific elastic strain energy  |-
        /*to be def.*/  double &sPD,                    // specific plastic dissipation    |---> Should be defined in Abaqus/Standard
        /*to be def.*/  double &sCD,                    // specific creep dissipation      |-
        //FOLLOWING ARGUMENTS: only in a fully coupled thermal-stress or thermal-electrical-structural analysis
        /*to be def.*/  double &rpl,                    // volumetric heat generation per unit time @ end of inc. caused by mech. working of material            
        /*to be def.*/  double ddSigma_ddTemp[],        // variation of stress with respect to temperature
        /*to be def.*/  double dRpl_dEpsilon[],         // variation of rpl with respect to strain increments
        /*to be def.*/  double &dRpl_dTemp,             // variation of rpl with respect to temperature
        const   double strain[],                        // array containing total strains @ beginning of increment     (only mechanical, no thermal); Shear Strain in engineering: e.g. gamma12 = 2*epsilon12
        const   double dStrain[],                       // array containing strain increment                           (only mechanical, no thermal)
        const   double time[2],                         // time[1]: value of step time @ beginning of current inc. or frequency
        //                                              // time[2]: value of total time @ beginning of current inc.
        const   double &dtime,                          // time increment
        const   double &temp,                           // temp @ start of inc.
        const   double &dTemp,                          // temp increment
        const   double preDef[],                        // array of interpolated values of predefined field variables @ this point @ start of inc., based on values read in nodes
        const   double dPreDef[],                       // array of inc. of pre. def. field variables
        const   char matName[80],                       // user defined material name, attention: If Intel Compiler is used, matNameLength *may*
                                                        // be passed directly after this argument
        const   int &nDirect,                           // number of direct stress components @ this point
        const   int &nShear,                            // number of engineering shear stress components @ this point
        const   int &nTensor,                           // size of stress and strain component array (nDirect + nShear)
        const   int &nStateVars,                            // number of solution dependent state variables associated with this mat. type
        const   double materialProperties[],                         // user def. array of mat. constants associated with this material
        const   int &nMaterialProperties,                            // number of user def. variables
        const   double coords[3],                       // coordinates of this point
        const   double dRot[9],                         // rotation increment matrix 3x3
        /*may be def.*/ double &pNewDT,                 // propagation for new time increment
        const   double &charElemLength,                 // characteristic element Length
        const   double dfGrd0[9],                       // deformation gradient @ beginning of increment      3x3 |
        const   double dfGrd1[9],                       // deformation gradient @ end of increment            3x3 |--> always stored as 3D-matrix
        const   int &noEl,                              // element number
        const   int &nPt,                               // integration Point number
        const   int &layer,                             // layer number (composite shells @ layered solids)
        const   int &kSectPt,                           // section point number within current layer
        const   int jStep[4],                           // step number **NOTE:** documentation and course material differ here, kStep[1] <-> jStep[4]
        //                                              // additional fields *may* be: procedure key, large deformation flag, perturbation step flag
        const   int &kInc,                              // increment Number
        const   int matNameLength                       // length of Material Name := 80, passed in when FORTRAN calls c/c++: Microsoft C compiler AND GCC (it *may* differ for IntelC++)
        ){       
          
        userLibrary::MaterialCode materialCode = static_cast<userLibrary::MaterialCode> ( stateVars[nStateVars-1] );
        if ( materialCode <= 0 ){
            const std::string materialName(matName);
            const std::string strippedName = materialName.substr(0, materialName.find_first_of(' ')). substr(0, materialName.find_first_of('-'));

            materialCode = userLibrary::getMaterialCodeFromName ( strippedName ); 

            stateVars[nStateVars-1] = static_cast<double> (materialCode);

        }

        BftMaterialHypoElastic* material = dynamic_cast<BftMaterialHypoElastic*> (bftMaterialFactory( materialCode, materialProperties, nMaterialProperties, noEl, nPt)); 

        const int nStateVarsForUmat = nStateVars - 1;

        if ( material->getNumberOfRequiredStateVars () > nStateVarsForUmat )  {
            const std::string materialName(matName);
            throw std::invalid_argument( MakeString() << "Material " << materialName.substr(0, materialName.find_first_of(' ')) << ": insufficient stateVars (" << nStateVars << " - 1) provided, but " << material->getNumberOfRequiredStateVars()<<" are required");
        }
        material->assignStateVars(stateVars, nStateVarsForUmat);

        material->setCharacteristicElementLength(charElemLength);
        
        double stress6[6], strain6[6], dStrain6[6], dStressDDStrain66[36];
        userLibrary::extendAbaqusToVoigt(stress6, stress, strain6, strain, dStrain6, dStrain, nDirect, nShear);

        if(nDirect == 3) 
            material->computeStress(stress6, dStressDDStrain66, strain6, dStrain6, time, dtime, pNewDT);
        else if(nDirect == 2)
            material->computePlaneStress(stress6, dStressDDStrain66, strain6, dStrain6, time, dtime, pNewDT);

        userLibrary::backToAbaqus(stress, stress6, dStressDDStrain, dStressDDStrain66, nDirect, nShear);

        delete material;
}

