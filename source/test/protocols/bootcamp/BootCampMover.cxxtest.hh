// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/bootcamp/BootCampMover.cxxtest.hh
/// @brief  A test suite for the BootCampMover, a Fold-Tree based MC mover made for the Rosetta Bootcamp.
/// @author decortja (joseph.a.decorte@vanderbilt.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>


// Project Headers


// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <string>

#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/moves/MoverFactory.hh>

static basic::Tracer TR("BootCampMover");


class BootCampMover : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}

    void test_MoverFactory() {
        protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover( "BootCampMover");
        protocols::bootcamp::BootCampMoverOP bcm_op = protocols::bootcamp::BootCampMoverOP( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( base_mover_op ) );
        TS_ASSERT_EQUALS( bcm_op->mover_name() , "BootCampMover");
    }

    void test_get_num_iterations() {
        protocols::bootcamp::BootCampMover test_bcm;
        test_bcm.set_num_iterations( 100);
        TS_ASSERT_EQUALS( test_bcm.get_num_iterations(), 100);
    }

    void test_get_sfxn() {
        protocols::bootcamp::BootCampMover test_bcm;
        core::scoring::ScoreFunctionOP new_sfxn_ = core::scoring::get_score_function();
        new_sfxn_->set_weight( core::scoring::overlap_chainbreak, 1);
        test_bcm.set_sfxn(  new_sfxn_);
        TS_ASSERT_EQUALS( test_bcm.get_sfxn(), new_sfxn_ );
    }

    // Parser testing


    void test_parse_niterations() {
        std::stringstream test_tag( "<BootCampMover niterations=100 scorefxn=\"testing123\"/>");
        utility::tag::TagCOP tag = utility::tag::Tag::create( test_tag);
        basic::datacache::DataMap test_datamap;
        core::scoring::ScoreFunctionOP test_sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction);
        test_datamap.add( "scorefxns", "testing123", test_sfxn);
        protocols::bootcamp::BootCampMover test_bcm;
        test_bcm.parse_my_tag( tag, test_datamap);
        TS_ASSERT_EQUALS(test_bcm.get_num_iterations(), 100);
    }

    void test_parse_sfxn() {
        std::stringstream test_tag( "<BootCampMover niterations=100 scorefxn=\"testing123\"/>");
        utility::tag::TagCOP tag = utility::tag::Tag::create( test_tag);
        basic::datacache::DataMap test_datamap;
        core::scoring::ScoreFunctionOP test_sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction);
        test_datamap.add( "scorefxns", "testing123", test_sfxn);
        protocols::bootcamp::BootCampMover test_bcm;
        test_bcm.parse_my_tag( tag, test_datamap);
        TS_ASSERT_EQUALS( test_bcm.get_sfxn(), test_sfxn);
    }
};