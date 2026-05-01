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

#include "AbaqusMarmotHelper.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

#include <Eigen/Dense>

// Instantiate the global tables for Abaqus/Standard exactly once in this translation unit
TableMap ElementTableMap;
TableMap MaterialTableMap;

extern "C" void uexternaldb_( int* LOP, int* LRESTART, double* TIME, double* DTIME, int* KSTEP, int* KINC )
{
  MutexInit( MutexID_UEL );
  
  if ( *LOP == 0 || *LOP == 4 ) {
    ElementTableMap.setLoaded( false );
    MaterialTableMap.setLoaded( false );
  }
}

// clang-format off
extern "C" void uel_(
  double rightHandSide[/*lVarx , nRightHandSide*/], // right hand side load vector(s) 1: common, 2: additional for RIKS (see documentation)
  double KMatrix[/*nDof * nDof*/],                  // stiffness matrix
  double stateVars[], // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
  double       energies[8], // may be updated: energies
  const int&   nDegreesOfFreedom,
  const int&   nRightHandSide,
  const int&   nStateVars,
  const double properties[/*nProperties*/],
  const int&   nProperties,
  const double coordinates[/*mcrd, nNodes*/], // undeformed coordinates of the node respective DOFs
  const int&   maxNCoords,                    // max number of coordinates (see documentation)
  const int&   nNodes,
  const double U[/*nDof*/],                             // current solution (end of increment)
  const double dU[/*mlvarx=nDof(?), nRightHandSide*/], // increment of solutions
  const double UDot[/*nDof*/],                          // first derivative (velocity..)
  const double UDotDot[/*nDof*/],                       // second derivative (acceleration..)
  const int&   elementType,                             // user defined element type id
  const double time[2],                                 // 1: time of step, 2: total time
  const double& dTime,                                  // time increment
  const int&   stepNumber,
  const int&   incrementNumber,
  const int&   elementNumber,
  const double solutionParams[],                   // solution procedure dependent parameters
  const int&   nDLoadActive,                       // id of load / flux currently active on this element
  const int distributedLoadTypes[/*mDLoads, * */], // An array containing the integers used to define distributed load types for the element.
  const double distributedLoadMags[/*mDloads, * */], // magnitudes @ end of increment
  const double Predef[/*nPredef*/],
  const int&   nPredef,
  const int    lFlags[],
  const int    mlvarx[],
  const double dDistributedLoadMags[/*mDloads, * */], // increment of magnitudes
  const int&   mDload, // total number of distributed loads and fluxes defined on this element
  double&      pNewDT,
  const int    integerProperties[],
  const int&   nIntegerProperties,
  const double& period )
// clang-format on
{
  try {
      loadIntToStringParameterTableOnceAndThreadSafe( "UEL_CODES", "UEL_ELEMENTS", ElementTableMap, MutexID_UEL );
      loadIntToStringParameterTableOnceAndThreadSafe( "UEL_CODES", "UEL_MATERIALS", MaterialTableMap, MutexID_UEL );

      const auto& elCodeToElName   = ElementTableMap;
      const auto& matCodeToMatName = MaterialTableMap;

      if ( nIntegerProperties < 3 ) {
        throw std::invalid_argument( std::format("Marmot: insufficient integer properties ({}) provided, at least 3 are required", nIntegerProperties) );
      }

      const int& elCode  = integerProperties[0];
      const int& matCode = integerProperties[1];
      const uint32_t additionalDefinitions = static_cast<uint32_t>(integerProperties[2]);

      const int nPropertiesMaterial = matCodeToMatName.at( matCode ).nProperties;
      const int nPropertiesElement  = nProperties - nPropertiesMaterial;

      const double* propertiesMaterial = &properties[0];
      const double* propertiesElement  = &properties[nPropertiesMaterial];

      auto theElement = std::unique_ptr< MarmotElement >(
        MarmotLibrary::MarmotElementFactory::createElement( elCodeToElName.at( elCode ).name, elementNumber ) );

      theElement->assignNodeCoordinates( coordinates );
      theElement->assignProperty( ElementProperties( propertiesElement, nPropertiesElement ) );
      theElement->assignProperty(
        MarmotMaterialSection( matCodeToMatName.at( matCode ).name, propertiesMaterial, nPropertiesMaterial ) );

      const int nNecessaryStateVars = theElement->getNumberOfRequiredStateVars();

      if ( nNecessaryStateVars > nStateVars ) {
        throw std::invalid_argument( std::format(
            "MarmotElement {} and material {}: insufficient stateVars ({}) provided, but {} are required",
            elCodeToElName.at(elCode).name, matCodeToMatName.at(matCode).name, nStateVars, nNecessaryStateVars) );
      }

      theElement->assignStateVars( stateVars, nStateVars );
      theElement->initializeYourself();

      int additionalDefinitionProperties = 0;
      
      if ( MainConstants::hasFlag(additionalDefinitions, MainConstants::AdditionalDefinitions::GeostaticStressDefinition) ) {
        if ( lFlags[0] == MainConstants::UelFlags1::GeostaticStress ) {
          const double* geostaticProperties = &propertiesElement[nPropertiesElement + additionalDefinitionProperties];
          theElement->setInitialConditions( MarmotElement::GeostaticStress, geostaticProperties );
        }
        additionalDefinitionProperties += 5;
      }

      if ( MainConstants::hasFlag(additionalDefinitions, MainConstants::AdditionalDefinitions::MarmotMaterialInitialization) && 
           stepNumber == 1 && incrementNumber == 1 ) {
        theElement->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
      }

      theElement->computeYourself( U, dU, rightHandSide, KMatrix, time, dTime, pNewDT );
      
  } catch (const std::exception& e) {
      handleAbaqusException(e, "UEL");
  } catch (...) {
      handleAbaqusUnknownException("UEL");
  }
}

// clang-format off
extern "C" void FOR_NAME(umat,UMAT)(
  /*to be def.*/ double stress[],    // stress vector in order: S11, S22, (S33), S12, (S13), (S23)
  /*to be def.*/ double stateVars[], // solution dependent state variables; passed in values @ beginning of increment -> set to values @ end of increment
  /*to be def.*/ double  dStress_dStrain[], // material Jacobian matrix ddSigma/ddEpsilon
  /*to be def.*/ double& sSE,               // specific elastic strain energy  |-
  /*to be def.*/ double& sPD, // specific plastic dissipation    |---> Should be defined in Abaqus/Standard
  /*to be def.*/ double& sCD, // specific creep dissipation      |-
  // FOLLOWING ARGUMENTS: only in a fully coupled thermal-stress or thermal-electrical-structural analysis
  /*to be def.*/ double& rpl, // volumetric heat generation per unit time @ end of inc. caused by mech. working of material
  /*to be def.*/ double  ddSigma_ddTemp[], // variation of stress with respect to temperature
  /*to be def.*/ double  dRpl_dEpsilon[],  // variation of rpl with respect to strain increments
  /*to be def.*/ double& dRpl_dTemp,       // variation of rpl with respect to temperature
  const double strain[],  // array containing total strains @ beginning of increment     (only mechanical, no thermal); Shear Strain in engineering: e.g. gamma12 = 2*epsilon12
  const double dStrain[], // array containing strain increment                           (only mechanical, no thermal)
  const double time[2],   // time[1]: value of step time @ beginning of current inc. or frequency // time[2]: value of total time @ beginning of current inc.
  const double& dtime,      // time increment
  const double& temp,       // temp @ start of inc.
  const double& dTemp,      // temp increment
  const double  preDef[],   // array of interpolated values of predefined field variables @ this point @ start of inc., based on values read in nodes
  const double dPreDef[],   // array of inc. of pre. def. field variables
  const char   matName[80], // user defined material name, attention: If Intel Compiler is used, matNameLength *may* be passed directly after this argument
  const int&             nDirect,    // number of direct stress components @ this point
  const int&             nShear,     // number of engineering shear stress components @ this point
  const int&             nTensor,    // size of stress and strain component array (nDirect + nShear)
  const int&             nStateVars, // number of solution dependent state variables associated with this mat. type
  const double           materialProperties[], // user def. array of mat. constants associated with this material
  const int&             nMaterialProperties,  // number of user def. variables
  const double           coords[3],            // coordinates of this point
  const double           dRot[9],              // rotation increment matrix 3x3
  /*may be def.*/ double& pNewDT,               // propagation for new time increment
  const double&           charElemLength,       // characteristic element Length
  const double           dfGrd0[9],            // deformation gradient @ beginning of increment     3x3 |
  const double dfGrd1[9], // deformation gradient @ end of increment            3x3 |--> always stored as 3D-matrix
  const int&   noEl,      // element number
  const int&   nPt,       // integration Point number
  const int&   layer,     // layer number (composite shells @ layered solids)
  const int&   kSectPt,   // section point number within current layer
  const int    jStep[4],  // step number **NOTE:** documentation and course material differ here, kStep[1] <-> jStep[4] // additional fields *may* be: procedure key, large deformation flag, perturbation step flag
  const int& kInc,        // increment Number
  const int matNameLength // length of Material Name := 80, passed in when FORTRAN calls c/c++: Microsoft C compiler AND
                          // GCC (it *may* differ for IntelC++)
)
// clang-format on
{
  try {
      std::string_view matNameView( matName, AbqStringLen );
      auto endPos = std::min(matNameView.find(' '), matNameView.find('-'));
      std::string strippedName(matNameView.substr(0, endPos));

      auto material = std::unique_ptr< MarmotMaterialHypoElastic >(
          MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( strippedName,
                                                                           materialProperties,
                                                                           nMaterialProperties,
                                                                           noEl ) );

      if ( material->getNumberOfRequiredStateVars() > nStateVars ) {
        throw std::invalid_argument( std::format(
            "MarmotMaterial {}: insufficient stateVars ({}) provided, but {} are required",
            strippedName, nStateVars, material->getNumberOfRequiredStateVars()) );
      }

      material->setCharacteristicElementLength( charElemLength );

      MarmotMaterialHypoElastic::timeInfo ti;
      ti.time = time[1];
      ti.dT   = dtime;  

      if ( nDirect == 3 ) {
        Eigen::Map<Eigen::VectorXd>       abqStress(stress, nTensor);
        Eigen::Map<const Eigen::VectorXd> abqDStrain(dStrain, nTensor);
        Eigen::Map<Eigen::MatrixXd>       abqDStressDStrain(dStress_dStrain, nTensor, nTensor);

        Eigen::Matrix<double, 6, 6> dStress_dStrain66 = Eigen::Matrix<double, 6, 6>::Zero();
        Eigen::Vector<double, 6>    dStrain6          = Eigen::Vector<double, 6>::Zero();

        // Stack-allocated array ensures zero heap allocations in the hot loop
        int abq2voigt[6] = {0}; 
        for ( int i = 0; i < nDirect; i++ ) abq2voigt[i] = i;
        for ( int i = 0; i < nShear; i++ )  abq2voigt[nDirect + i] = 3 + i;

        MarmotMaterialHypoElastic::state3D state;
        state.stateVars           = stateVars;
        state.strainEnergyDensity = sSE;
        state.stress.setZero(); 

        for ( int i = 0; i < nTensor; i++ ) {
          state.stress(abq2voigt[i]) = abqStress(i);
          dStrain6(abq2voigt[i])     = abqDStrain(i);
        }

        material->computeStress( state, dStress_dStrain66.data(), dStrain6.data(), ti );

        sSE = state.strainEnergyDensity;

        for ( int i = 0; i < nTensor; i++ ) {
          abqStress(i) = state.stress(abq2voigt[i]);
          for ( int j = 0; j < nTensor; j++ ) {
            abqDStressDStrain(i, j) = dStress_dStrain66(abq2voigt[i], abq2voigt[j]);
          }
        }
      }
      else if ( nDirect == 2 ) {
        MarmotMaterialHypoElastic::state2D state;
        state.stateVars           = stateVars;
        state.strainEnergyDensity = sSE;
        
        for (int i = 0; i < nTensor; i++) {
            state.stress(i) = stress[i];
        }

        material->computePlaneStress( state, dStress_dStrain, dStrain, ti );

        sSE = state.strainEnergyDensity;
        
        for (int i = 0; i < nTensor; i++) {
            stress[i] = state.stress(i);
        }
      }
      else if ( nDirect == 1 ) {
        MarmotMaterialHypoElastic::state1D state;
        state.stateVars           = stateVars;
        state.strainEnergyDensity = sSE;
        
        state.stress = stress[0]; 

        material->computeUniaxialStress( state, dStress_dStrain, dStrain, ti );

        sSE       = state.strainEnergyDensity;
        stress[0] = state.stress;
      }
      
  } catch (const std::exception& e) {
      handleAbaqusException(e, "UMAT");
  } catch (...) {
      handleAbaqusUnknownException("UMAT");
  }
}
