/* ---------------------------------------------------------------------
 *                                       _   
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_ 
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_ 
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 * 
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck, 
 * 2020 - today
 * 
 * festigkeitslehre@uibk.ac.at
 * 
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 * 
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#include <aba_for_c.h>
#include "Marmot/Marmot.h"
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMaterialHypoElastic.h"

namespace MainConstants
{
    enum AdditionalDefinitions {
        GeostaticStressDefiniton = 0x01 << 0,
        MarmotMaterialInitialization = 0x01 << 1,
    };

    enum UelFlags1 {         
        GeostaticStress = 61, // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to                               // procedure types.     
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
        void FOR_NAME(stdb_abqerr, STDB_ABQERR)(const int *lop, const char* stringZT, const int *intArray, const double *realArray, const char *appendix, const int lengthString, const int lengthAppendix);
        void FOR_NAME(xit, XIT)();
}


extern "C" void FOR_NAME(uel, UEL)(
        double rightHandSide[/*lVarx , nRightHandSide*/],           // right hand side load vector(s) 1: common, 2: additional for RIKS (see documentation)
        double KMatrix[/*nDof * nDof*/],                            // stiffness matrix 
        double stateVars[],                                         // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        double energies[8],                                         // may be updated: energies
        const int &nDegreesOfFreedom ,
        const int &nRightHandSide,                
        const int &nStateVars, 
        const double properties[/*nProperties*/],
        const int &nProperties,
        const double coordinates[/*mcrd, nNodes*/],                 // undeformed coordinates of the node respective DOFs
        const int &maxNCoords,                                      // max number of coordinates (see documentation)
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
            throw std::invalid_argument( MakeString() << "Marmot: insufficient integer properties (" << nIntegerProperties <<") provided, but 5 are required");

        MarmotLibrary::ElementCode elementCode =  static_cast<MarmotLibrary::ElementCode> ( integerProperties[0]);
        MarmotLibrary::MaterialCode materialID =  static_cast<MarmotLibrary::MaterialCode>( integerProperties[1] );
        const int nPropertiesElement =          integerProperties[2];
        const int nPropertiesUmat =             integerProperties[3];
        const int additionalDefinitions =       integerProperties[4];

        const double* propertiesUmat =    &properties[0];
        const double* propertiesElement = &properties[nPropertiesUmat];

        auto theElement = std::unique_ptr<MarmotElement> ( MarmotLibrary::MarmotElementFactory::createElement(elementCode,  elementNumber) );

        theElement->assignProperty( ElementProperties( propertiesElement, nPropertiesElement ) );

        theElement->assignProperty( MarmotMaterialSection( materialID, propertiesUmat, nPropertiesUmat) );

        const int nNecessaryStateVars = theElement->getNumberOfRequiredStateVars();

        if ( nNecessaryStateVars > nStateVars )
            throw std::invalid_argument( MakeString() << "MarmotElement with code " << elementCode << " and material " << materialID << ": insufficient stateVars (" << nStateVars <<") provided, but "<< nNecessaryStateVars <<" are required");

        theElement->assignStateVars(stateVars, nStateVars);

        theElement->initializeYourself(coordinates);

        int additionalDefinitionProperties = 0;
        if( additionalDefinitions & MainConstants::AdditionalDefinitions::GeostaticStressDefiniton ) {
            if(lFlags[0] == MainConstants::UelFlags1::GeostaticStress){
                const double* geostaticProperties = &propertiesElement[ nPropertiesElement + additionalDefinitionProperties ];
                theElement->setInitialConditions( MarmotElement::GeostaticStress, geostaticProperties ); }
            additionalDefinitionProperties += 5;  
        }

        if( additionalDefinitions & MainConstants::AdditionalDefinitions::MarmotMaterialInitialization 
                && stepNumber == 1 && incrementNumber == 1) {
            theElement->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr); 
        }

        // compute K and P 
        theElement->computeYourself(U , dU, rightHandSide, KMatrix, time, dTime, pNewDT); 

        //// compute distributed loads in nodal forces and add it to P 
        //for (int i =0; i<mDload; i++){
            //if ([i]<1.e-16)
                //continue;
            //theElement->computeDistributedLoad(MarmotElement::Pressure, rightHandSide, distributedLoadTypes[i], &distributedLoadMags[i], time, dTime);}

}

extern "C" void FOR_NAME(umat, UMAT)(
        /*to be def.*/  double stress[],                // stress vector in order: S11, S22, (S33), S12, (S13), (S23) 
        /*to be def.*/  double stateVars[],             // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
        /*to be def.*/  double dStressDDStrain[],       // material Jacobian matrix ddSigma/ddEpsilon
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
        const   int &nStateVars,                        // number of solution dependent state variables associated with this mat. type
        const   double materialProperties[],            // user def. array of mat. constants associated with this material
        const   int &nMaterialProperties,               // number of user def. variables
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
          
        MarmotLibrary::MaterialCode materialCode = static_cast<MarmotLibrary::MaterialCode> ( stateVars[nStateVars-1] );
        if ( materialCode <= 0 ){
            const std::string materialName(matName);
            const std::string strippedName = materialName.substr(0, materialName.find_first_of(' ')). substr(0, materialName.find_first_of('-'));

            materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName ( strippedName ); 

            stateVars[nStateVars-1] = static_cast<double> (materialCode);
        }

        auto material = std::unique_ptr<MarmotMaterialHypoElastic> (
                dynamic_cast<MarmotMaterialHypoElastic*> (
                    MarmotLibrary::MarmotMaterialFactory::createMaterial( materialCode, materialProperties, nMaterialProperties, noEl )));

        const int nStateVarsForUmat = nStateVars - 1;

        if ( material->getNumberOfRequiredStateVars () > nStateVarsForUmat )  {
            const std::string materialName(matName);
            throw std::invalid_argument( MakeString() << "MarmotMaterial " << materialName.substr(0, materialName.find_first_of(' ')) 
                    << ": insufficient stateVars (" << nStateVars << " - 1) provided, but " << material->getNumberOfRequiredStateVars()<<" are required");
        }
        material->assignStateVars(stateVars, nStateVarsForUmat);

        material->setCharacteristicElementLength(charElemLength);
        
        double stress6[6] = {};
        double dStrain6[6] = {};
        double dStressDDStrain66[36] = {};

        int abq2voigt [nTensor];
        for (int i = 0; i < nDirect; i++)
            abq2voigt[ i ] = i;
        for (int i = 0; i < nShear; i++)
            abq2voigt[ nDirect+i ]   = 3 + i;
        
        // expand Voigt
        for (int i = 0; i < nTensor; i++) {
            stress6 [ abq2voigt[ i ] ] = stress[ i ];
            dStrain6 [ abq2voigt[ i ] ] = dStrain [ i ];
        }

	    // call material
        if(nDirect == 3) 
            material->computeStress(stress6, dStressDDStrain66,  dStrain6, time, dtime, pNewDT);
        else if(nDirect == 2)
            material->computePlaneStress(stress6, dStressDDStrain66,  dStrain6, time, dtime, pNewDT);

	    // condense Voigt
        for (int i = 0; i < nTensor; i++) {
            stress[ i ] = stress6[ abq2voigt[ i ] ];
            for (int j = 0; j < nTensor; j++)
                dStressDDStrain[ nTensor * i +  j ] = dStressDDStrain66[ 6 * abq2voigt[ i ] + abq2voigt [ j ] ];
            }
}
