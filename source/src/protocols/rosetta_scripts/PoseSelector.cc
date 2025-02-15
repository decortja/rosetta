// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/rosetta_scripts/PoseSelector.hh
/// @brief  Base class for pose selectors used by MultiplePoseMover
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_rosetta_scripts_PoseSelector_CC
#define INCLUDED_protocols_rosetta_scripts_PoseSelector_CC

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelector.hh>

// Package headers

// Project headers

// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers

static basic::Tracer TR( "protocols.rosetta_scripts.PoseSelector" );

namespace protocols {
namespace rosetta_scripts {

PoseSelector::PoseSelector() = default;

PoseSelector::~PoseSelector() = default;

void PoseSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &
)
{
	TR << "***WARNING!!!! WARNING!!!*** parse_my_tag has been invoked for this PoseSelector but it hasn't been defined. Are you sure this is appropriate?" << std::endl;
}

} // rosetta_scripts
} // protocols

#endif //INCLUDED_protocols_rosetta_scripts_PoseSelector_CC
