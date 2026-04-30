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
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include <aba_for_c.h>
#include <SMAAspUserSubroutines.h>

#include <Eigen/Dense>

#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <format>
#include <vector>
#include <ranges>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>

namespace MainConstants {
  enum AdditionalDefinitions {
    GeostaticStressDefiniton     = 0x01 << 0,
    MarmotMaterialInitialization = 0x01 << 1,
  };

  enum UelFlags1 {
    GeostaticStress = 61, // Geostatic stress field according to Abaqus Analysis User's Guide Tab. 5.1.2-1 Keys to procedure types.
  };
} // namespace MainConstants

enum MutexIDs {
  MutexID_UEL = 1,
};

extern "C" {
// clang-format off
void FOR_NAME(stdb_abqerr,STDB_ABQERR)
  // clang-format on
  ( const int*    lop,
    const char*   stringZT,
    const int*    intArray,
    const double* realArray,
    const char*   appendix,
    const int     lengthString,
    const int     lengthAppendix );
// clang-format off
void FOR_NAME(xit,XIT)();
// clang-format on
}

extern "C" {
void settablecollection_( char* tcName, int* jError, int tcName_len );

void queryparametertable_( char* parameterTableLabel,
                           int*  numParams,
                           int*  numRows,
                           int*  jError,
                           int   parameterTableLabel_len );

void getparametertablerow_( char*   parameterTableLabel,
                            int*    jRow,
                            int*    numParams,
                            int*    iParamsDataType,
                            int*    iParams,
                            double* rParams,
                            char*   cParams,
                            int*    jError,
                            int     parameterTableLabel_len,
                            int     cParams_len );
}

// C++23 String helper for 80-char Fortran buffers
void make_fstr80(char* dest, std::string_view src) {
    std::ranges::fill(dest, dest + 80, ' ');
    std::ranges::copy(src.substr(0, std::min<size_t>(80, src.size())), dest);
}

// Modernized TableMap using composition to avoid std container inheritance issues
struct MarmotInfo {
  std::string name;
  int         nProperties;
};

class TableMap {
public:
  [[nodiscard]] bool isLoaded() const noexcept { return loaded; }
  void setLoaded(bool val) noexcept { loaded = val; }

  MarmotInfo& operator[](int key) { return data[key]; }
  const MarmotInfo& at(int key) const { return data.at(key); }
  void clear() { data.clear(); }

private:
  bool loaded = false;
  std::unordered_map<int, MarmotInfo> data;
};

TableMap ElementTableMap;
TableMap MaterialTableMap;

// RAII Wrapper for Abaqus Mutex to guarantee exception safety
class AbaqusScopedLock {
public:
    explicit AbaqusScopedLock(int mutexId) : id(mutexId) { MutexLock(id); }
    ~AbaqusScopedLock() { MutexUnlock(id); }
    
    // Prevent copying
    AbaqusScopedLock(const AbaqusScopedLock&) = delete;
    AbaqusScopedLock& operator=(const AbaqusScopedLock&) = delete;
private:
    int id;
};

void readAllMarmotInfoInto( TableMap& map, std::string_view tableCollectionName, std::string_view tableLabel )
{
  char tcName[80], tblName[80];
  make_fstr80( tcName, tableCollectionName );
  make_fstr80( tblName, tableLabel );

  int jError = 0;

  // ---------------------------------------------------
  // Activate table collection
  // ---------------------------------------------------
  settablecollection_( tcName, &jError, 80 );
  if ( jError != 0 ) {
    std::cout << std::format("ERROR: Cannot activate table collection {}\n", tableCollectionName);
    return;
  }

  // ---------------------------------------------------
  // Query number of rows and parameters
  // ---------------------------------------------------
  int numParams = 0, numRows = 0;
  queryparametertable_( tblName, &numParams, &numRows, &jError, 80 );

  if ( jError != 0 || numRows == 0 ) {
    std::cout << std::format("ERROR: Cannot query parameter table {}\n", tableLabel);
    return;
  }

  // ---------------------------------------------------
  // Prepare buffers safely with std::vector (RAII)
  // ---------------------------------------------------
  const int maxParams   = numParams;
  const int cParams_len = 80 * maxParams;

  std::vector<int>    iParamsDataType(maxParams);
  std::vector<int>    iParams(maxParams);
  std::vector<double> rParams(maxParams);
  std::string         cParams(cParams_len, ' '); // initializes directly with spaces

  map.clear();

  std::cout << std::format("Reading parameter table {} with {} rows and {} parameters per row.\n", 
                           tableLabel, numRows, numParams);

  // ---------------------------------------------------
  // Loop over all rows
  // ---------------------------------------------------
  for ( int row = 1; row <= numRows; row++ ) {
    getparametertablerow_( tblName,
                           &row,
                           &numParams,
                           iParamsDataType.data(),
                           iParams.data(),
                           rParams.data(),
                           cParams.data(),
                           &jError,
                           80,
                           cParams_len );

    if ( jError != 0 ) {
      std::cout << std::format("ERROR: getParameterTableRow failed at row {}\n", row);
      continue;
    }

    // integer is param 1
    int key         = iParams[0];
    int nProperties = iParams[2];

    // string is param 2 → index 1 → starts at offset 80*(index). Zero-allocation trim.
    std::string_view strView(cParams.data() + 80, 80);
    auto lastChar = strView.find_last_not_of(' ');
    std::string str = (lastChar == std::string_view::npos) ? "" : std::string(strView.substr(0, lastChar + 1));

    map[key] = MarmotInfo{ str, nProperties };

    std::cout << std::format("Loaded row {}: {} -> {}, number of properties = {}\n", 
                             row, key, str, nProperties);
  }
}

void loadIntToStringParameterTableOnceAndThreadSafe( std::string_view collection,
                                                     std::string_view label,
                                                     TableMap&        map )
{
  if ( map.isLoaded() )
    return;

  AbaqusScopedLock lock(MutexID_UEL);

  if ( map.isLoaded() ) {
    return;
  }
  
  readAllMarmotInfoInto( map, collection, label );
  map.setLoaded( true );
}

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
  loadIntToStringParameterTableOnceAndThreadSafe( "UEL_CODES", "UEL_ELEMENTS", ElementTableMap );
  loadIntToStringParameterTableOnceAndThreadSafe( "UEL_CODES", "UEL_MATERIALS", MaterialTableMap );

  const auto& elCodeToElName   = ElementTableMap;
  const auto& matCodeToMatName = MaterialTableMap;

  if ( nIntegerProperties < 3 )
    throw std::invalid_argument( std::format("Marmot: insufficient integer properties ({}) provided, at least 3 are required", nIntegerProperties) );

  const int& elCode                = integerProperties[0];
  const int& matCode               = integerProperties[1];
  const int& additionalDefinitions = integerProperties[2];

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

  if ( nNecessaryStateVars > nStateVars )
    throw std::invalid_argument( std::format(
        "MarmotElement {} and material {}: insufficient stateVars ({}) provided, but {} are required",
        elCodeToElName.at(elCode).name, matCodeToMatName.at(matCode).name, nStateVars, nNecessaryStateVars) );

  theElement->assignStateVars( stateVars, nStateVars );
  theElement->initializeYourself();

  int additionalDefinitionProperties = 0;
  if ( additionalDefinitions & MainConstants::AdditionalDefinitions::GeostaticStressDefiniton ) {
    if ( lFlags[0] == MainConstants::UelFlags1::GeostaticStress ) {
      const double* geostaticProperties = &propertiesElement[nPropertiesElement + additionalDefinitionProperties];
      theElement->setInitialConditions( MarmotElement::GeostaticStress, geostaticProperties );
    }
    additionalDefinitionProperties += 5;
  }

  if ( additionalDefinitions & MainConstants::AdditionalDefinitions::MarmotMaterialInitialization && stepNumber == 1 &&
       incrementNumber == 1 ) {
    theElement->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  theElement->computeYourself( U, dU, rightHandSide, KMatrix, time, dTime, pNewDT );
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
  // Zero-allocation parsing
  std::string_view matNameView( matName, 80 );
  auto endPos = std::min(matNameView.find(' '), matNameView.find('-'));
  std::string strippedName(matNameView.substr(0, endPos));

  auto material = std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( strippedName,
                                                                       materialProperties,
                                                                       nMaterialProperties,
                                                                       noEl ) );

  const int nStateVarsForUmat = nStateVars;

  if ( material->getNumberOfRequiredStateVars() > nStateVarsForUmat ) {
    throw std::invalid_argument( std::format(
        "MarmotMaterial {}: insufficient stateVars ({}) provided, but {} are required",
        strippedName, nStateVars, material->getNumberOfRequiredStateVars()) );
  }

  material->setCharacteristicElementLength( charElemLength );

  // Set up the new timeInfo struct
  MarmotMaterialHypoElastic::timeInfo ti;
  ti.time = time[1]; // Total time at beginning of increment
  ti.dT   = dtime;   // Time increment

  if ( nDirect == 3 ) {
    // Map Abaqus raw pointers directly to Eigen types (Zero-copy, Column-Major)
    Eigen::Map<Eigen::VectorXd>       abqStress(stress, nTensor);
    Eigen::Map<const Eigen::VectorXd> abqDStrain(dStrain, nTensor);
    Eigen::Map<Eigen::MatrixXd>       abqDStressDStrain(dStress_dStrain, nTensor, nTensor);

    // Initialize 6x6 and 6x1 Eigen containers for the material library
    Eigen::Matrix<double, 6, 6> dStress_dStrain66 = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Vector<double, 6>    dStrain6          = Eigen::Vector<double, 6>::Zero();

    std::vector<int> abq2voigt(nTensor);
    for ( int i = 0; i < nDirect; i++ ) abq2voigt[i] = i;
    for ( int i = 0; i < nShear; i++ )  abq2voigt[nDirect + i] = 3 + i;

    // Set up the 3D state struct
    MarmotMaterialHypoElastic::state3D state;
    state.stateVars           = stateVars;
    state.strainEnergyDensity = sSE;
    state.stress.setZero(); 

    // Expand Voigt notation
    for ( int i = 0; i < nTensor; i++ ) {
      state.stress(abq2voigt[i]) = abqStress(i);
      dStrain6(abq2voigt[i])     = abqDStrain(i);
    }

    // Call the updated compute function
    material->computeStress( state, dStress_dStrain66.data(), dStrain6.data(), ti );

    // Update specific elastic strain energy
    sSE = state.strainEnergyDensity;

    // Condense Voigt notation
    for ( int i = 0; i < nTensor; i++ ) {
      abqStress(i) = state.stress(abq2voigt[i]);
      for ( int j = 0; j < nTensor; j++ ) {
        abqDStressDStrain(i, j) = dStress_dStrain66(abq2voigt[i], abq2voigt[j]);
      }
    }
  }
  else if ( nDirect == 2 ) {
    Eigen::Map<Eigen::MatrixXd> abqDStressDStrain(dStress_dStrain, nTensor, nTensor);

    // Set up the 2D state struct
    MarmotMaterialHypoElastic::state2D state;
    state.stateVars           = stateVars;
    state.strainEnergyDensity = sSE;
    
    // Safest way to populate a fixed-size vector from a runtime-sized pointer
    for (int i = 0; i < nTensor; i++) {
        state.stress(i) = stress[i];
    }

    material->computePlaneStress( state, abqDStressDStrain.data(), dStrain, ti );

    sSE = state.strainEnergyDensity;
    
    // Write back to Abaqus array
    for (int i = 0; i < nTensor; i++) {
        stress[i] = state.stress(i);
    }
  }
  else if ( nDirect == 1 ) {
    Eigen::Map<Eigen::MatrixXd> abqDStressDStrain(dStress_dStrain, 1, 1);

    // Set up the 1D state struct
    MarmotMaterialHypoElastic::state1D state;
    state.stateVars           = stateVars;
    state.strainEnergyDensity = sSE;
    
    state.stress = stress[0]; 

    material->computeUniaxialStress( state, abqDStressDStrain.data(), dStrain, ti );

    sSE       = state.strainEnergyDensity;
    stress[0] = state.stress;
  }
}
