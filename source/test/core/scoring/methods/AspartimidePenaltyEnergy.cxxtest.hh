// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/AspartimidePenaltyEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::AspartimidePenaltyEnergy
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers


// Package Headers
#include <test/core/init_util.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/annotated_sequence.hh>

#include <basic/Tracer.hh>



#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/DeclareBond.hh>

//Auto Headers



#include <core/pose/Pose.hh> // AUTO IWYU For Pose
#include <core/scoring/Energies.hh> // AUTO IWYU For Energies
#include <core/scoring/ScoreFunction.hh> // AUTO IWYU For ScoreFunction


static basic::Tracer TR("core.scoring.methods.AspartimidePenaltyEnergy.cxxtest");

// --------------- Test Class --------------- //


class AspartimidePenaltyEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_cyclic_geometry_scoring() {
		TR << "Staring AspartimidePenaltyEnergyTests::test_cyclic_geometry_scoring()." << std::endl;

		// Set up test pose and cyclize
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");
		protocols::simple_moves::DeclareBond connect_ends;
		connect_ends.set( 8, "C", 1, "N", false );
		connect_ends.apply(pose);

		// Set up score function
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::aspartimide_penalty, 10.0 );

		// Score with no aspartate:
		(sfxn)(pose);
		TR << "\tSEQ\tEXP\tACT" << std::endl;
		TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);

		{ //Scope 1: one aspartate
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(8);
			mutres.set_res_name("ASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 2: two aspartates, one with DASP-LALA sequence.
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(6);
			mutres.set_res_name("DASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 3: two aspartates, one aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("DASN");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 4: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(1);
			mutres.set_res_name("SER");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 5: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(1);
			mutres.set_res_name("THR");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 6: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(1);
			mutres.set_res_name("GLY");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 7: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(8);
			mutres.set_res_name("DASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 8: two aspartates, one aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(1);
			mutres.set_res_name("DILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 9: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(1);
			mutres.set_res_name("ILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

	}

	void test_L_scoring() {
		TR << "Staring AspartimidePenaltyEnergyTests::test_L_scoring()." << std::endl;

		// Set up test pose
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");

		// Set up score function
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::aspartimide_penalty, 10.0 );

		// Score with no aspartate:
		(sfxn)(pose);
		TR << "\tSEQ\tEXP\tACT" << std::endl;
		TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);

		{ //Scope 1: one aspartate
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(4);
			mutres.set_res_name("ASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 2: two aspartates
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(6);
			mutres.set_res_name("ASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 3: two aspartates, one aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("ASN");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 4: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(5);
			mutres.set_res_name("SER");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 5: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("THR");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 6: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(5);
			mutres.set_res_name("GLY");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 7: one aspartate, one aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(4);
			mutres.set_res_name("ILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 8: one aspartate, no aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("ILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 9: one aspartate, one aspartimide sequence with D-residue
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("DILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
			TR << "\tResidue 7 is " << pose.residue(7).name() << ".  The is_d_aa() function returns " << (pose.residue(7).type().is_d_aa() ? "true." : "false.") << std::endl;
		}

	}


	void test_D_scoring() {
		TR << "Staring AspartimidePenaltyEnergyTests::test_D_scoring()." << std::endl;

		// Set up test pose
		core::pose::Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");
		for ( core::Size i=1, imax=pose.size(); i<=imax; ++i ) {
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(i);
			mutres.set_res_name("DALA");
			mutres.apply(pose);
		}

		// Set up score function
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::aspartimide_penalty, 10.0 );

		// Score with no aspartate:
		(sfxn)(pose);
		TR << "\tSEQ\tEXP\tACT" << std::endl;
		TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);

		{ //Scope 1: one aspartate
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(4);
			mutres.set_res_name("DASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 2: two aspartates
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(6);
			mutres.set_res_name("DASP");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 3: two aspartates, one aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("DASN");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 4: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(5);
			mutres.set_res_name("DSER");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 5: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("DTHR");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 6: two aspartates, two aspartimide sequences
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(5);
			mutres.set_res_name("GLY");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 500 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 500.0, 1e-6);
		}

		{ //Scope 7: one aspartate, one aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(4);
			mutres.set_res_name("DILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
		}

		{ //Scope 8: one aspartate, no aspartimide sequence
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("DILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 0 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 0.0, 1e-6);
		}

		{ //Scope 9: one aspartate, one aspartimide sequence with D-residue
			protocols::simple_moves::MutateResidue mutres;
			mutres.set_target(7);
			mutres.set_res_name("ILE");
			mutres.apply(pose);
			(sfxn)(pose);
			TR << "\t" << pose.sequence() << "\t" << 250 << "\t" << pose.energies().total_energy() << std::endl;
			TS_ASSERT_DELTA(pose.energies().total_energy(), 250.0, 1e-6);
			TR << "\tResidue 7 is " << pose.residue(7).name() << ".  The is_d_aa() function returns " << (pose.residue(7).type().is_d_aa() ? "true." : "false.") << std::endl;
		}

	}

};


