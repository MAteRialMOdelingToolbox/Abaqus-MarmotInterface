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
#pragma once

#include <SMAAspUserSubroutines.h>
#include <aba_for_c.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <format>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Global constant for Fortran string buffers
constexpr size_t AbqStringLen = 80;

namespace MainConstants {
  // Strongly typed 32-bit bitmask for initial conditions and extra definitions
  enum class AdditionalDefinitions : uint32_t {
    None                         = 0,
    GeostaticStressDefinition    = 1u << 0,
    MarmotMaterialInitialization = 1u << 1,
    UnusedFlag02                 = 1u << 2,
    UnusedFlag03                 = 1u << 3,
    UnusedFlag04                 = 1u << 4,
    UnusedFlag05                 = 1u << 5,
    UnusedFlag06                 = 1u << 6,
    UnusedFlag07                 = 1u << 7,
    UnusedFlag08                 = 1u << 8,
    UnusedFlag09                 = 1u << 9,
    UnusedFlag10                 = 1u << 10,
    UnusedFlag11                 = 1u << 11,
    UnusedFlag12                 = 1u << 12,
    UnusedFlag13                 = 1u << 13,
    UnusedFlag14                 = 1u << 14,
    UnusedFlag15                 = 1u << 15,
    UnusedFlag16                 = 1u << 16,
    UnusedFlag17                 = 1u << 17,
    UnusedFlag18                 = 1u << 18,
    UnusedFlag19                 = 1u << 19,
    UnusedFlag20                 = 1u << 20,
    UnusedFlag21                 = 1u << 21,
    UnusedFlag22                 = 1u << 22,
    UnusedFlag23                 = 1u << 23,
    UnusedFlag24                 = 1u << 24,
    UnusedFlag25                 = 1u << 25,
    UnusedFlag26                 = 1u << 26,
    UnusedFlag27                 = 1u << 27,
    UnusedFlag28                 = 1u << 28,
    UnusedFlag29                 = 1u << 29,
    UnusedFlag30                 = 1u << 30,
    UnusedFlag31                 = 1u << 31
  };

  // Helper for clean bitwise checking
  [[nodiscard]] constexpr bool hasFlag( uint32_t bitmask, AdditionalDefinitions flag ) noexcept
  {
    return ( bitmask & static_cast< uint32_t >( flag ) ) != 0;
  }

  enum UelFlags1 {
    GeostaticStress = 61,
  };
} // namespace MainConstants

enum MutexIDs {
  MutexID_UEL   = 1,
  MutexID_VUMAT = 2, // Added for future Explicit expansion
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

// ALL functions below must be marked 'inline' to prevent multiple definition linker errors

// C++23 String helper for 80-char Fortran buffers
inline void make_fstr80( char* dest, std::string_view src )
{
  std::ranges::fill( dest, dest + AbqStringLen, ' ' );
  std::ranges::copy( src.substr( 0, std::min< size_t >( AbqStringLen, src.size() ) ), dest );
}

/**
 * Universal Abaqus Output Router
 * LOP =  1: Informational message
 * LOP = -1: Warning message
 * LOP = -2: Error message (continue)
 * LOP = -3: Error message (stop immediately)
 */
inline void printAbaqusMessage( std::string_view msg, int lop )
{
  int    dummyInt  = 0;
  double dummyReal = 0.0;

  // Truncate to Abaqus's 500 character limit to prevent buffer overflows
  std::string msgStr = std::string( msg.substr( 0, 500 ) );

  FOR_NAME( stdb_abqerr, STDB_ABQERR )( &lop, msgStr.c_str(), &dummyInt, &dummyReal, "", msgStr.length(), 0 );
}

inline void handleAbaqusException( const std::exception& e, std::string_view routineName )
{
  std::string msg = std::format( "MARMOT FATAL ERROR in {}: {}", routineName, e.what() );
  printAbaqusMessage( msg, -3 ); // -3 = Stop immediately
}

inline void handleAbaqusUnknownException( std::string_view routineName )
{
  std::string msg = std::format( "MARMOT FATAL ERROR: Unknown exception caught in {}.", routineName );
  printAbaqusMessage( msg, -3 ); // -3 = Stop immediately
}

// Modernized TableMap using composition
struct MarmotInfo {
  std::string name;
  int         nProperties;
};

class TableMap {
public:
  [[nodiscard]] bool isLoaded() const noexcept { return loaded; }
  void               setLoaded( bool val ) noexcept { loaded = val; }

  MarmotInfo&       operator[]( int key ) { return data[key]; }
  const MarmotInfo& at( int key ) const { return data.at( key ); }
  void              clear() { data.clear(); }

private:
  bool                                  loaded = false;
  std::unordered_map< int, MarmotInfo > data;
};

// RAII Wrapper for Abaqus Mutex to guarantee exception safety
class AbaqusScopedLock {
public:
  explicit AbaqusScopedLock( int mutexId ) : id( mutexId ) { MutexLock( id ); }
  ~AbaqusScopedLock() { MutexUnlock( id ); }

  // Prevent copying
  AbaqusScopedLock( const AbaqusScopedLock& )            = delete;
  AbaqusScopedLock& operator=( const AbaqusScopedLock& ) = delete;

private:
  int id;
};

inline void readAllMarmotInfoInto( TableMap& map, std::string_view tableCollectionName, std::string_view tableLabel )
{
  char tcName[80], tblName[80];
  make_fstr80( tcName, tableCollectionName );
  make_fstr80( tblName, tableLabel );

  int jError = 0;
  settablecollection_( tcName, &jError, AbqStringLen );
  if ( jError != 0 ) {
    printAbaqusMessage( std::format( "MARMOT WARNING: Cannot activate table collection {}", tableCollectionName ), -1 );
    return;
  }

  int numParams = 0, numRows = 0;
  queryparametertable_( tblName, &numParams, &numRows, &jError, AbqStringLen );
  if ( jError != 0 || numRows == 0 ) {
    printAbaqusMessage( std::format( "MARMOT WARNING: Cannot query parameter table {}", tableLabel ), -1 );
    return;
  }

  const int maxParams   = numParams;
  const int cParams_len = AbqStringLen * maxParams;

  std::vector< int >    iParamsDataType( maxParams );
  std::vector< int >    iParams( maxParams );
  std::vector< double > rParams( maxParams );
  std::string           cParams( cParams_len, ' ' );

  map.clear();

  printAbaqusMessage( std::format( "MARMOT INFO: Reading parameter table {} with {} rows and {} parameters per row.",
                                   tableLabel,
                                   numRows,
                                   numParams ),
                      1 );

  for ( int row = 1; row <= numRows; row++ ) {
    getparametertablerow_( tblName,
                           &row,
                           &numParams,
                           iParamsDataType.data(),
                           iParams.data(),
                           rParams.data(),
                           cParams.data(),
                           &jError,
                           AbqStringLen,
                           cParams_len );

    if ( jError != 0 ) {
      printAbaqusMessage( std::format( "MARMOT WARNING: getParameterTableRow failed at row {}", row ), -1 );
      continue;
    }

    int key         = iParams[0];
    int nProperties = iParams[2];

    std::string_view strView( cParams.data() + AbqStringLen, AbqStringLen );
    auto             lastChar = strView.find_last_not_of( ' ' );
    std::string str = ( lastChar == std::string_view::npos ) ? "" : std::string( strView.substr( 0, lastChar + 1 ) );

    map[key] = MarmotInfo{ str, nProperties };

    printAbaqusMessage( std::format( "MARMOT INFO: Loaded row {}: {} -> {}, number of properties = {}",
                                     row,
                                     key,
                                     str,
                                     nProperties ),
                        1 );
  }
}

inline void loadIntToStringParameterTableOnceAndThreadSafe( std::string_view collection,
                                                            std::string_view label,
                                                            TableMap&        map,
                                                            int              mutexId )
{
  if ( map.isLoaded() )
    return;

  AbaqusScopedLock lock( mutexId );

  if ( map.isLoaded() )
    return;

  readAllMarmotInfoInto( map, collection, label );
  map.setLoaded( true );
}
