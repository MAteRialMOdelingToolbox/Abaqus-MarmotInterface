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
// clang-format off
#ifndef NO_ABAQUS
#include <aba_for_c.h>
#endif
#ifndef FOR_NAME
#define FOR_NAME(a, b) a##_
#endif
// clang-format on
#include "Marmot/MarmotExplicitLibrary.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticExplicit.h"
#include <Eigen/Core>
#include <memory>
#include <string>

// clang-format off
extern "C" void FOR_NAME(vumat, VUMAT)(
// clang-format on
  const int&    nBlocks,
  const int&    nDirect,
  const int&    nShear,
  const int&    nStateVars,
  const int&    nfieldv,
  const int&    nMaterialProperties,
  const int&    lanneal,
  const double& stepTime,
  const double& totalTime,
  double*       dTArray,
  const char*   matName,
  const double* coordMp,
  const double* charLength,
  const double* materialProperties,
  const double* density,
  const double* dStrain,
  const double* relSpinInc,
  const double* tempOld,
  const double* stretchOld,
  const double* defgradOld,
  const double* fieldOld,
  const double* stressOld,
  const double* stateVarsOld,
  const double* enerInternOld,
  const double* enerInelasOld,
  const double* tempNew,
  const double* stretchNew,
  const double* defgradNew,
  const double* fieldNew,
  double*       stressNew,
  double*       stateVarsNew,
  double*       enerInternNew,
  double*       enerInelasNew,
  const int     matNameLength // length of Material Name := 80, passed in when FORTRAN calls c/c++:
                              // Microsoft C compiler AND GCC (it *may* differ for IntelC++)
)

{

  using namespace Eigen;

  const std::string materialName( matName );
  const std::string strippedName = materialName.substr( 0, materialName.find_first_of( ' ' ) )
                                     .substr( 0, materialName.find_first_of( '-' ) );

  int materialCode = MarmotLibrary::MarmotMaterialExplicitFactory::getMaterialCodeFromName( strippedName );

  auto material = std::unique_ptr< MarmotMaterialHypoElasticExplicit >(
    dynamic_cast< MarmotMaterialHypoElasticExplicit* >(
      MarmotLibrary::MarmotMaterialExplicitFactory::createMaterial( materialCode,
                                                                    materialProperties,
                                                                    nMaterialProperties,
                                                                    0 ) ) );

  // alternative way to create the material (doesn't add a lot of performance)
  /* auto material = Marmot::Materials::CDPExplicit( materialProperties, nMaterialProperties, 0 ); */

  const int nTensor = nDirect + nShear;

  // Map to stress old and new, and initialize new with old
  Map< const MatrixXd > stressOldBlock( stressOld, nBlocks, nTensor );
  Map< MatrixXd >       stressNewBlock( stressNew, nBlocks, nTensor );
  stressNewBlock = stressOldBlock;

  // Map to dStrain, and multiply shear terms by a factor of two
  MatrixXd dStrainBlock = Map< const MatrixXd >( dStrain, nBlocks, nTensor );

  // we need to multiply all shear terms by a factor of two in Abaqus explicit:
  dStrainBlock.block( 0, nDirect, nBlocks, nShear ) *= 2.0;

  // Map to stateVars old and new, and initialize new with old
  Map< const MatrixXd > stateVarsOldBlock( stateVarsOld, nBlocks, nStateVars );
  Map< MatrixXd >       stateVarsNewBlock( stateVarsNew, nBlocks, nStateVars );
  stateVarsNewBlock = stateVarsOldBlock;

  const double    dTime = dTArray[0];
  Map< VectorXd > shearModuliForDT( dTArray + 1, nBlocks );
  Map< VectorXd > bulkModuliForDT( dTArray + 1 + nBlocks, nBlocks );

  Map< const VectorXd > materialProperties_( materialProperties, nMaterialProperties );

  Map< VectorXd >       internalEnergyDensityBlock( enerInternNew, nBlocks );
  Map< VectorXd >       dissipatedEnergyDensityBlock( enerInelasNew, nBlocks );
  Map< const VectorXd > densityBlock( density, nBlocks );
  Map< const VectorXd > charLengthBlock( charLength, nBlocks );

  if ( nShear == 3 ) {
    // Furthermore, we need to swap the last two shear components in Abaqus explicit (only in 3d):
    // dStrain
    VectorXd temp         = dStrainBlock.col( 4 );
    dStrainBlock.col( 4 ) = dStrainBlock.col( 5 );
    dStrainBlock.col( 5 ) = temp;

    // stress
    temp                    = stressNewBlock.col( 4 );
    stressNewBlock.col( 4 ) = stressNewBlock.col( 5 );
    stressNewBlock.col( 5 ) = temp;

    // same for stress:

    MatrixXd stateVarsNewBlock_RowMajor = stateVarsNewBlock.transpose();
    MatrixXd dStrainBlock3d_RowMajor    = dStrainBlock.transpose();
    MatrixXd stressBlock3D_RowMajor     = stressNewBlock.transpose();

    for ( int b = 0; b < nBlocks; b++ ) {

      /* material->assignStateVars( stateVarsNewBlock_RowMajor.col(b).data(), nStateVars ); */

      material->computeStress( stressBlock3D_RowMajor.col( b ).data(),
                               stateVarsNewBlock_RowMajor.col( b ).data(),
                               nStateVars,
                               internalEnergyDensityBlock.data() + b,
                               dissipatedEnergyDensityBlock.data() + b,
                               dStrainBlock3d_RowMajor.col( b ).data(),
                               densityBlock( b ),
                               totalTime,
                               dTime,
                               shearModuliForDT( b ),
                               bulkModuliForDT( b ),
                               charLengthBlock( b ) );
    }

    stateVarsNewBlock = stateVarsNewBlock_RowMajor.transpose();
    stressNewBlock    = stressBlock3D_RowMajor.transpose();

    // we need to swap back the last two shear components in Abaqus explicit (only in 3d):
    temp                    = stressNewBlock.col( 4 );
    stressNewBlock.col( 4 ) = stressNewBlock.col( 5 );
    stressNewBlock.col( 5 ) = temp;
  }

  else if ( nShear == 1 ) {

    // plane strain case

    // make a full 3d dStrainBlock:
    MatrixXd dStrainBlock3d( nBlocks, 6 );
    dStrainBlock3d.setZero();

    dStrainBlock3d.block( 0, 0, nBlocks, 4 ) = dStrainBlock;

    // same for stress:
    MatrixXd stressBlock3D( nBlocks, 6 );
    stressBlock3D.setZero();
    stressBlock3D.block( 0, 0, nBlocks, 4 ) = stressNewBlock;

    MatrixXd stateVarsNewBlock_RowMajor = stateVarsNewBlock.transpose();
    MatrixXd dStrainBlock3d_RowMajor    = dStrainBlock3d.transpose();
    MatrixXd stressBlock3D_RowMajor     = stressBlock3D.transpose();

    for ( int b = 0; b < nBlocks; b++ ) {

      /* material->assignStateVars( stateVarsNewBlock_RowMajor.col(b).data(), nStateVars ); */

      material->computeStress( stressBlock3D_RowMajor.col( b ).data(),
                               stateVarsNewBlock_RowMajor.col( b ).data(),
                               nStateVars,
                               internalEnergyDensityBlock.data() + b,
                               dissipatedEnergyDensityBlock.data() + b,
                               dStrainBlock3d_RowMajor.col( b ).data(),
                               densityBlock( b ),
                               totalTime,
                               dTime,
                               shearModuliForDT( b ),
                               bulkModuliForDT( b ),
                               charLengthBlock( b ) );
    }

    stateVarsNewBlock = stateVarsNewBlock_RowMajor.transpose();

    // condense 3D stress to plane strain stress:
    stressNewBlock = stressBlock3D_RowMajor.transpose().block( 0, 0, nBlocks, 4 );
  }
}
// clang-format off
extern "C" void FOR_NAME(vuel, VUEL)
// clang-format on
                                      ( const int&    nBlock,
                                        double*       rhs,
                                        double*       amass,
                                        double*       dTStable,
                                        double*       stateVars,
                                        const int&    nStateVars,
                                        double*       energies,
                                        const int&    nNodes,
                                        const int&    nDofElement,
                                        const double* properties,
                                        const int&    nProperties,
                                        const int*    integerProperties,
                                        const int&    nIntegerProperties,
                                        const double* coordinates,
                                        const int&    mcrd,
                                        const double* U,
                                        const double* dU,
                                        const double* UDot,
                                        const double* UDotDot,
                                        const int&    jtype,
                                        const int*    jElem,
                                        const double* time,
                                        const double* period,
                                        const double& dTime,
                                        const double& dTimePrev,
                                        const int&    kstep,
                                        const int&    kinc,
                                        const int*    lflags,
                                        const double* massScaleFactor, // this parameter seems to be always zero; needs to be
                                                                       // checked
                                        const double* predef,
                                        const int&    npredef,
                                        const int&    jdltyp,
                                        const double* adlmag )
{
  if ( nIntegerProperties != 5 )
    throw std::invalid_argument( MakeString()
                                 << "Marmot: insufficient integer properties (" << nIntegerProperties
                                 << ") provided, but 5 are required: elementCode, materialCode, nPropertiesElement, "
                                    "nPropertiesMaterial, additionalDefinitions, nStateVarsMaterial " );

  const int elementCode           = integerProperties[0];
  const int materialCode          = integerProperties[1];
  const int nPropertiesMaterial   = integerProperties[2];
  const int nPropertiesElement    = integerProperties[3];
  const int additionalDefinitions = integerProperties[4];

  const double* materialProperties = properties;
  const double* elementProperties  = properties + nPropertiesMaterial;

  const int nElEnergies = 12;

  using namespace Eigen;

  Map< MatrixXd > stateVarsBlock( stateVars, nBlock, nStateVars );
  Map< MatrixXd > rhsBlock( rhs, nBlock, nDofElement );
  Map< MatrixXd > amassBlock( amass, nBlock, nDofElement * nDofElement );
  Map< MatrixXd > energiesBlock( energies, nBlock, nElEnergies );

  Map< const MatrixXd > UBlock( U, nBlock, nDofElement );
  Map< const MatrixXd > dUBlock( dU, nBlock, nDofElement );
  Map< const MatrixXd > UDotBlock( UDot, nBlock, nDofElement );
  Map< const MatrixXd > UDotDotBlock( UDotDot, nBlock, nDofElement );
  Map< const MatrixXd > coordinatesBlock( coordinates, nBlock, mcrd * nNodes );

  MatrixXd stateVarsBlock_RowMajor = stateVarsBlock.transpose();

  auto theMaterial = ( MarmotLibrary::MarmotMaterialExplicitFactory::createMaterial( materialCode,
                                                                                     materialProperties,
                                                                                     nPropertiesMaterial,
                                                                                     0 ) );

  auto theElement = std::unique_ptr< MarmotElementExplicit >(
    MarmotLibrary::MarmotElementExplicitFactory::createElement( elementCode, 0 ) );

  theElement->assignProperties( elementProperties, nPropertiesElement );

  theElement->assignMaterial( theMaterial );

  const auto& procedureType = lflags[0]; // 17 = explicit dynamic, 74 = explicit fully coupled thermo-mechanical
  const auto& nlgeom        = lflags[1]; // 0 = linear, 1 = nonlinear
  const auto& opCode        = lflags[2]; // 1 = mass matrix
                                         // 2 = internal forces and critical time step
                                         // 3 = external forces
  const double currentTime = time[1];    // 0 = step time, 1 = total time

  const MatrixXd coordinatesBlock_RowMajor = coordinatesBlock.transpose();
  amassBlock.setZero();

  MatrixXd elCoordinates_RowMajor( mcrd, nNodes );

  if ( opCode == 1 ) {
    MatrixXd amassBlock_RowMajor = amassBlock.transpose();

    for ( int b = 0; b < nBlock; b++ ) {

      // Unlike Abaqus/Standard, Abaqus/Explicit stores the node coordinates in colmajor (nNode,mcrd) (instead of
      // colmajor (mcrd,nNode)) so we have to swap again here
      elCoordinates_RowMajor = Map< const MatrixXd >( coordinatesBlock_RowMajor.col( b ).data(), nNodes, mcrd )
                                 .transpose();

      theElement->assignNodeCoordinates( elCoordinates_RowMajor.data() );
      theElement->initializeYourself();
      theElement->computeConsistentMassMatrix( amassBlock_RowMajor.col( b ).data() );
      theElement->lumpMassMatrix( amassBlock_RowMajor.col( b ).data() );
    }

    amassBlock = amassBlock_RowMajor.transpose();
  }

  else if ( opCode == 2 ) {

    const MatrixXd UBlock_RowMajor       = UBlock.transpose();
    const MatrixXd dUBlock_RowMajor      = dUBlock.transpose();
    const MatrixXd UDotBlock_RowMajor    = UDotBlock.transpose();
    const MatrixXd UDotDotBlock_RowMajor = UDotDotBlock.transpose();

    MatrixXd rhsBlock_RowMajor       = rhsBlock.transpose();
    MatrixXd energiesBlock_RowMajor  = energiesBlock.transpose();
    MatrixXd stateVarsBlock_RowMajor = stateVarsBlock.transpose();

    for ( int b = 0; b < nBlock; b++ ) {

      elCoordinates_RowMajor = Map< const MatrixXd >( coordinatesBlock_RowMajor.col( b ).data(), nNodes, mcrd )
                                 .transpose();
      theElement->assignNodeCoordinates( elCoordinates_RowMajor.data() );
      theElement->initializeYourself();

      theElement->computeKernels(
        // input fields
        UBlock_RowMajor.col( b ).data(),
        dUBlock_RowMajor.col( b ).data(),
        UDotBlock_RowMajor.col( b ).data(),
        UDotDotBlock_RowMajor.col( b ).data(),

        // output
        rhsBlock_RowMajor.col( b ).data(),
        energiesBlock_RowMajor.col( b ).data(),
        stateVarsBlock_RowMajor.col( b ).data(),

        nStateVars,

        currentTime,
        dTime,
        dTStable[b] );
    }

    rhsBlock       = rhsBlock_RowMajor.transpose();
    energiesBlock  = energiesBlock_RowMajor.transpose();
    stateVarsBlock = stateVarsBlock_RowMajor.transpose();
  }

  delete theMaterial;
}
