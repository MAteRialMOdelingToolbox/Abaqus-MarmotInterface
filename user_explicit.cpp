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
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialHypoElasticExplicit.h"
#include <Eigen/Core>
#include <aba_for_c.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

extern "C" void FOR_NAME( vumat, VUMAT )(
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
  const int matNameLength // length of Material Name := 80, passed in when FORTRAN calls c/c++: Microsoft C compiler AND
                          // GCC (it *may* differ for IntelC++)
)

{
  const std::string materialName( matName );
  const std::string strippedName = materialName.substr( 0, materialName.find_first_of( ' ' ) )
                                     .substr( 0, materialName.find_first_of( '-' ) );

  int materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( strippedName );

  auto material = std::unique_ptr< MarmotMaterialHypoElasticExplicit >(
    dynamic_cast< MarmotMaterialHypoElasticExplicit* >(
      MarmotLibrary::MarmotMaterialFactory::createMaterial( materialCode,
                                                            materialProperties,
                                                            nMaterialProperties,
                                                            0 ) ) );

  // alternative way to create the material (doesn't add a lot of performance)
  /* auto material = Marmot::Materials::CDPExplicit( materialProperties, nMaterialProperties, 0 ); */

  const int nTensor = nDirect + nShear;

  // Map to stress old and new, and initialize new with old
  Eigen::Map< const Eigen::MatrixXd > stressOldBlock( stressOld, nBlocks, nTensor );
  Eigen::Map< Eigen::MatrixXd >       stressNewBlock( stressNew, nBlocks, nTensor );
  stressNewBlock = stressOldBlock;

  // Map to dStrain, and multiply shear terms by a factor of two
  Eigen::MatrixXd dStrainBlock = Eigen::Map< const Eigen::MatrixXd >( dStrain, nBlocks, nTensor );

  // we need to multiply all shear terms by a factor of two in Abaqus explicit:
  dStrainBlock.block( 0, nDirect, nBlocks, nShear ) *= 2.0;

  // Map to stateVars old and new, and initialize new with old
  Eigen::Map< const Eigen::MatrixXd > stateVarsOldBlock( stateVarsOld, nBlocks, nStateVars );
  Eigen::Map< Eigen::MatrixXd >       stateVarsNewBlock( stateVarsNew, nBlocks, nStateVars );
  stateVarsNewBlock = stateVarsOldBlock;

  const double                  dTime = dTArray[0];
  Eigen::Map< Eigen::VectorXd > shearModuliForDT( dTArray + 1, nBlocks );
  Eigen::Map< Eigen::VectorXd > bulkModuliForDT( dTArray + 1 + nBlocks, nBlocks );

  Eigen::Map< const Eigen::VectorXd > materialProperties_( materialProperties, nMaterialProperties );

  Eigen::Map< Eigen::VectorXd >       internalEnergyDensityBlock( enerInternNew, nBlocks );
  Eigen::Map< Eigen::VectorXd >       dissipatedEnergyDensityBlock( enerInelasNew, nBlocks );
  Eigen::Map< const Eigen::VectorXd > densityBlock( density, nBlocks );
  Eigen::Map< const Eigen::VectorXd > charLengthBlock( charLength, nBlocks );

  if ( nShear == 3 ) {
    // Furthermore, we need to swap the last two shear components in Abaqus explicit (only in 3d):
    // dStrain
    Eigen::VectorXd temp  = dStrainBlock.col( 4 );
    dStrainBlock.col( 4 ) = dStrainBlock.col( 5 );
    dStrainBlock.col( 5 ) = temp;

    // stress
    temp                    = stressNewBlock.col( 4 );
    stressNewBlock.col( 4 ) = stressNewBlock.col( 5 );
    stressNewBlock.col( 5 ) = temp;

    // same for stress:

    Eigen::MatrixXd stateVarsNewBlock_RowMajor = stateVarsNewBlock.transpose();
    Eigen::MatrixXd dStrainBlock3d_RowMajor    = dStrainBlock.transpose();
    Eigen::MatrixXd stressBlock3D_RowMajor     = stressNewBlock.transpose();

    for ( int b = 0; b < nBlocks; b++ ) {

      /* material->assignStateVars( stateVarsNewBlock_RowMajor.col(b).data(), nStateVars ); */

      material->computeStress( stressBlock3D_RowMajor.col( b ).data(),
                               stateVarsNewBlock_RowMajor.col( b ).data(),
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
    Eigen::MatrixXd dStrainBlock3d( nBlocks, 6 );
    dStrainBlock3d.setZero();

    dStrainBlock3d.block( 0, 0, nBlocks, 4 ) = dStrainBlock;

    // same for stress:
    Eigen::MatrixXd stressBlock3D( nBlocks, 6 );
    stressBlock3D.setZero();
    stressBlock3D.block( 0, 0, nBlocks, 4 ) = stressNewBlock;

    Eigen::MatrixXd stateVarsNewBlock_RowMajor = stateVarsNewBlock.transpose();
    Eigen::MatrixXd dStrainBlock3d_RowMajor    = dStrainBlock3d.transpose();
    Eigen::MatrixXd stressBlock3D_RowMajor     = stressBlock3D.transpose();

    for ( int b = 0; b < nBlocks; b++ ) {

      /* material->assignStateVars( stateVarsNewBlock_RowMajor.col(b).data(), nStateVars ); */

      material->computeStress( stressBlock3D_RowMajor.col( b ).data(),
                               stateVarsNewBlock_RowMajor.col( b ).data(),
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
