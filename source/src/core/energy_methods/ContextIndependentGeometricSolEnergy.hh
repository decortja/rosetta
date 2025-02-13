// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.hh
/// @brief  Geometric solvation energy.
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_ContextIndependentGeometricSolEnergy_hh

// Unit Headers
#include <core/energy_methods/ContextIndependentGeometricSolEnergy.fwd.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh>
#include <core/types.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>


namespace core {
namespace energy_methods {


class ContextIndependentGeometricSolEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;
public:


	ContextIndependentGeometricSolEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	/// @brief copy c-tor
	ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// attempt to precalculate backbone/backbone energies in advance
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const & ) const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & scfxn ) const override;

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		bool res_moving_wrt_eachother
	) const override;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const res1,
		Size const res2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &
	) const;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;

	virtual
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &
	) const;

	bool
	use_extended_residue_pair_energy_interface() const override;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData & pair_data
	) const override;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & ires,
		conformation::Residue const & jres,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const &,
		pose::Pose const & pose,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
	) const override;

	bool
	requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const override;

	void
	setup_for_minimizing_for_residue(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		basic::datacache::BasicDataCache &,
		core::scoring::ResSingleMinimizationData & res_data_cache
	) const override;

	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	/// This evaluates everything for now,
	/// but eventually may want to split this
	/// based on backbone/backbone vs. others,
	/// as is carried out in HBondEnergy.cc
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const override;

	// Undefined, commenting out to fix PyRosetta build
	/* void
	eval_atom_derivative_intra_RNA(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;
	*/


	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	/// f1 and f2 are zeroed
	// virtual
	// void
	// eval_atom_derivative(
	//  id::AtomID const & atom_id,
	//  pose::Pose const & pose,
	//  kinematics::DomainMap const &,
	//  core::scoring::ScoreFunction const &,
	//  core::scoring::EnergyMap const & weights,
	//  Vector & F1,
	//  Vector & F2
	// ) const;


	bool
	defines_intrares_energy( core::scoring::EnergyMap const & weights ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & ,
		core::scoring::EnergyMap & emap
	) const override;

	Distance
	atomic_interaction_cutoff() const override;

	core::Size version() const override;

	/// @brief GeometricSolEnergy is context sensitive
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const override;

	void
	precalculate_bb_bb_energy_for_design(
		pose::Pose const &
	) const;


private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	core::scoring::methods::EnergyMethodOptions const & options_;

	scoring::geometric_solvation::GeometricSolEnergyEvaluatorOP evaluator_;

	mutable Real precalculated_bb_bb_energy_;

	mutable bool using_extended_method_;


};

} // scoring
} // core

#endif
