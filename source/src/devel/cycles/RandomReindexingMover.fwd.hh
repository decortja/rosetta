// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_devel_cycles_RandomReindexingMover_fwd_hh
#define INCLUDED_devel_cycles_RandomReindexingMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace cycles {

class RandomReindexingMover;
typedef utility::pointer::owning_ptr <RandomReindexingMover> RandomReindexingMoverOP;
typedef utility::pointer::owning_ptr <RandomReindexingMover const> RandomReindexingMoverCOP;

}
}

#endif
