// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   platform/windows/msvc/platform/types.hh
/// @brief  Platform-specific types for Windows using MSVC
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_platform_windows_msvc_platform_types_hh
#define INCLUDED_platform_windows_msvc_platform_types_hh


// Windows SDK headers
#include <basetsd.h>

// C++ headers
#include <stdint.h>
#include <cstddef>

/// @brief Fixed size types
// int64_t; // 64-bit signed integer
// uint64_t; // 64-bit unsigned integer


/// @brief Scalable size types
// intptr_t  Pointer-sized signed integer
// uintptr_t  Pointer-sized unsigned integer
// ssize_t  Signed size


namespace platform {

typedef std::size_t  Size;
typedef SSIZE_T      SSize;
typedef std::size_t  uint;

// Floating point precision control scalar
#ifdef ROSETTA_FLOAT // Real == float
typedef  float  Real;
#else // Real == double
typedef  double  Real;
#endif

namespace file {


/// @brief Are file/path names case sensitive?
bool const CASE_SENSITIVE( false );

/// @brief Volume specifier used?
bool const VOLUME_USED( true );

/// @brief File path separator
char const PATH_SEPARATOR( '\\' );


} // namespace file


} // namespace platform


#endif // INCLUDED_platform_types_HH
