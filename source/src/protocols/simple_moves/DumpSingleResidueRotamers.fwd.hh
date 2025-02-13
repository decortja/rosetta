// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DumpSingleResidueRotamers.fwd.hh
/// @brief Given a residue index, dump all of the rotamers to individual PDB files within 0-1 sd of the mean
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_DumpSingleResidueRotamers_fwd_hh
#define INCLUDED_protocols_simple_moves_DumpSingleResidueRotamers_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace simple_moves {

class DumpSingleResidueRotamers;

typedef utility::pointer::shared_ptr< DumpSingleResidueRotamers > DumpSingleResidueRotamersOP;
typedef utility::pointer::shared_ptr< DumpSingleResidueRotamers const > DumpSingleResidueRotamersCOP;

} //protocols
} //simple_moves

#endif //INCLUDED_protocols_simple_moves_DumpSingleResidueRotamers_fwd_hh
