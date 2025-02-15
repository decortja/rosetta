// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/DirectReadoutEnergy.hh
/// @brief  Statistically derived DNA contact potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_energy_methods_DirectReadoutEnergy_hh
#define INCLUDED_core_energy_methods_DirectReadoutEnergy_hh

// Unit Headers
#include <core/energy_methods/DirectReadoutEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>


// Utility headers

namespace core {
namespace energy_methods {


/// @brief  Implementation of Kono and Sarai's knowledge-based protein-DNA interaction energy
/// @details  Could be a CI2B, but interaction threshold is large, so in the short term defining as
/// WholeStructure energy.

class DirectReadoutEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;
public:


	DirectReadoutEnergy();


	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////


	/// @brief  Implementation which is currently not used
	///
	void
	my_residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const;

	using parent::finalize_total_energy;

	/// @brief  All the work happens here
	virtual
	void
	finalize_total_energy(
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::scoring::EnergyMap & emap
	) const;


	/// @brief  No graphs required.
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {};


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	core::scoring::dna::DirectReadoutPotential const & potential_;
	core::Size version() const override;

};


}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
