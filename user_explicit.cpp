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
// clang-format off
#ifndef NO_ABAQUS
#include <aba_for_c.h>
#endif
#ifndef FOR_NAME
#define FOR_NAME(a, b) a##_
#endif
// clang-format on
#define ABAQUS_EXPLICIT
#include "AbaqusMarmotHelper.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include <Eigen/Core>
#include <memory>
#include <string>

// Instantiate the global tables for Abaqus/Standard exactly once in this translation unit
TableMap ElementTableMap;
TableMap MaterialTableMap;

extern "C" void FOR_NAME( vexternaldb,
                          VEXTERNALDB )( int* lOp, int* i_Array, int* niArray, double* r_Array, int* nrArray )
{
  // --- 0-Based Indexing Conversions (Fortran is 1-based) ---
  // Contents of i_Array
  const int i_int_nTotalNodes      = 0;
  const int i_int_nTotalElements   = 1;
  const int i_int_kStep            = 2;
  const int i_int_kInc             = 3;
  const int i_int_iStatus          = 4;
  const int i_int_lWriteRestart    = 5;
  const int i_int_ExtraOutputFrame = 6;

  // Possible values for the lOp argument
  const int j_int_StartAnalysis  = 0;
  const int j_int_StartStep      = 1;
  const int j_int_SetupIncrement = 2;
  const int j_int_StartIncrement = 3;
  const int j_int_EndIncrement   = 4;
  const int j_int_EndStep        = 5;
  const int j_int_EndAnalysis    = 6;

  // Possible values for i_Array[i_int_iStatus]
  const int j_int_Continue          = 0;
  const int j_int_TerminateStep     = 1;
  const int j_int_TerminateAnalysis = 2;

  // Contents of r_Array
  const int i_flt_TotalTime = 0;
  const int i_flt_StepTime  = 1;
  const int i_flt_dTime     = 2;

  // --- Accessing variables (Dereferencing Pointers) ---
  int operation = *lOp;
  int kStep     = i_Array[i_int_kStep];
  int kInc      = i_Array[i_int_kInc];

  if ( operation == j_int_StartAnalysis ) {

    ElementTableMap.setLoaded( false );
    MaterialTableMap.setLoaded( false );
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

  try {
    loadIntToStringParameterTableOnceAndThreadSafe( "VUEL_CODES", "VUEL_ELEMENTS", ElementTableMap );
    loadIntToStringParameterTableOnceAndThreadSafe( "VUEL_CODES", "VUEL_MATERIALS", MaterialTableMap );

    const auto& elCodeToElName   = ElementTableMap;
    const auto& matCodeToMatName = MaterialTableMap;

    if ( nIntegerProperties < 3 ) {
      throw std::invalid_argument(
        std::format( "Marmot: insufficient integer properties ({}) provided, at least 3 are required",
                     nIntegerProperties ) );
    }

    const int&     elCode                = integerProperties[0];
    const int&     matCode               = integerProperties[1];
    const uint32_t additionalDefinitions = static_cast< uint32_t >( integerProperties[2] );

    const int nPropertiesMaterial = matCodeToMatName.at( matCode ).nProperties;
    const int nPropertiesElement  = nProperties - nPropertiesMaterial;

    const double* propertiesMaterial = &properties[0];
    const double* propertiesElement  = &properties[nPropertiesMaterial];

    auto theElement = std::unique_ptr< MarmotElement >(
      MarmotLibrary::MarmotElementFactory::createElement( elCodeToElName.at( elCode ).name, jElem[0] ) );

    theElement->assignProperty( ElementProperties( propertiesElement, nPropertiesElement ) );
    theElement->assignProperty(
      MarmotMaterialSection( matCodeToMatName.at( matCode ).name, propertiesMaterial, nPropertiesMaterial ) );

    const int nNecessaryStateVars = theElement->getNumberOfRequiredStateVars();

    if ( nNecessaryStateVars > nStateVars ) {
      throw std::invalid_argument(
        std::format( "MarmotElement {} and material {}: insufficient stateVars ({}) provided, but {} are required",
                     elCodeToElName.at( elCode ).name,
                     matCodeToMatName.at( matCode ).name,
                     nStateVars,
                     nNecessaryStateVars ) );
    }

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

    const auto& procedureType = lflags[0];
    const auto& nlgeom        = lflags[1];
    const auto& opCode        = lflags[2];

    const MatrixXd coordinatesBlock_RowMajor = coordinatesBlock.transpose();
    amassBlock.setZero();

    MatrixXd elCoordinates_RowMajor( mcrd, nNodes );

    if ( opCode == 1 ) {
      MatrixXd amassBlock_RowMajor = amassBlock.transpose();

      for ( int b = 0; b < nBlock; b++ ) {

        elCoordinates_RowMajor = Map< const MatrixXd >( coordinatesBlock_RowMajor.col( b ).data(), nNodes, mcrd )
                                   .transpose();

        theElement->assignNodeCoordinates( elCoordinates_RowMajor.data() );
        theElement->assignStateVars( stateVarsBlock_RowMajor.col( b ).data(), nStateVars );
        theElement->initializeYourself();

        // Adapted to use the new inertia methods. Explicit dynamics typically requires the lumped mass.
        theElement->computeConsistentInertia( amassBlock_RowMajor.col( b ).data() );

        // lump the matrix:
        MatrixXd massMatrix = Map< MatrixXd >( amassBlock_RowMajor.col( b ).data(), nDofElement, nDofElement );
        VectorXd lumpedMass = massMatrix.rowwise().sum();
        massMatrix.setZero();
        massMatrix.diagonal() = lumpedMass;

        amassBlock_RowMajor.col( b ) = Map< VectorXd >( massMatrix.data(), nDofElement * nDofElement );
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

        // State variables must be assigned explicitly now before computation
        theElement->assignStateVars( stateVarsBlock_RowMajor.col( b ).data(), nStateVars );
        theElement->initializeYourself();

        // dTStable requires an isolated variable since the signature uses double&
        double stableTimeStep = dTStable[b];

        // Replaces computeKernels. UDot and UDotDot are omitted as they are not present in the new explicit interface.
        theElement->computeYourselfExplicit( UBlock_RowMajor.col( b ).data(),
                                             dUBlock_RowMajor.col( b ).data(),
                                             rhsBlock_RowMajor.col( b ).data(),
                                             time,
                                             dTime,
                                             stableTimeStep );

        dTStable[b] = stableTimeStep;

        // Extract internal energy. Abaqus expects ALLIE (internal energy) at index 0 of the energies array.
        double internalEnergy = 0.0;
        theElement->computeInternalEnergy( internalEnergy );

        /* C     energy array indices */
        /*       parameter ( iElPd = 1, */
        /*      *            iElCd = 2, */
        /*      *            iElIe = 3, */
        /*      *            iElTs = 4, */
        /*      *            iElDd = 5, */
        /*      *            iElBv = 6, */
        /*      *            iElDe = 7, */
        /*      *            iElHe = 8, */
        /*      *            iUnused = 9, */
        /*      *            iElTh = 10, */
        /*      *            iElDmd = 11, */
        /*      *            iElDc = 12, */
        /*      *            nElEnergy = 12) */

        energiesBlock_RowMajor.col( b )( 2 ) = internalEnergy;
      }

      rhsBlock       = -rhsBlock_RowMajor.transpose();
      energiesBlock  = energiesBlock_RowMajor.transpose();
      stateVarsBlock = stateVarsBlock_RowMajor.transpose();
    }
  }
  catch ( const std::exception& e ) {
    handleAbaqusException( e, "VUEL" );
  }
  catch ( ... ) {
    handleAbaqusUnknownException( "VUEL" );
  }
}
