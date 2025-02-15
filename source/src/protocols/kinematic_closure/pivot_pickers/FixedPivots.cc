// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedPivots.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

using namespace std;

FixedPivots::FixedPivots(core::Size start, core::Size stop, core::Size cut) {
	start_ = start;
	stop_ = stop;
	cut_ = cut;
}

Loop FixedPivots::pick(Pose const &, Loop const &) {
	return Loop(start_, stop_, cut_);
}

}
}
}
