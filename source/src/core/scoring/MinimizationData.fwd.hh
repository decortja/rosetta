// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/MinimizationData.fwd.hh
/// @brief  Forward declaration of the container class for use by certain EnergyMethods
//          during derivative- and score-function evaluation within minimization routines.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_MinimizationData_fwd_hh
#define INCLUDED_core_scoring_MinimizationData_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class ResSingleMinimizationData;
class ResPairMinimizationData;

typedef utility::pointer::shared_ptr< ResSingleMinimizationData > ResSingleMinimizationDataOP;
typedef utility::pointer::shared_ptr< ResSingleMinimizationData const > ResSingleMinimizationDataCOP;

typedef utility::pointer::shared_ptr< ResPairMinimizationData > ResPairMinimizationDataOP;
typedef utility::pointer::shared_ptr< ResPairMinimizationData const > ResPairMinimizationDataCOP;

}
}

#endif
