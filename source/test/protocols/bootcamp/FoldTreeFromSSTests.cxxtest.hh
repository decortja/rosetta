// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/FoldTreeFromSS.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/vector1.hh>

/// Project headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <test/util/pose_funcs.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// C++ headers
#include <string>

//Auto Headers


// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

//    protocols::bootcamp::identify_secondary_structure_spans( std::string const & ss_string )
//    {
//        utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
//        core::Size strand_start = -1;
//        for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
//            if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
//                if ( int( strand_start ) == -1 ) {
//                    strand_start = ii;
//                } else if ( ss_string[ii] != ss_string[strand_start] ) {
//                    ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
//                    strand_start = ii;
//                }
//            } else {
//                if ( int( strand_start ) != -1 ) {
//                    ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
//                    strand_start = -1;
//                }
//            }
//        }
//        if ( int( strand_start ) != -1 ) {
//            // last residue was part of a ss-element
//            ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
//        }
//        for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
//            std::cout << "SS Element " << ii << " from residue "
//                      << ss_boundaries[ ii ].first << " to "
//                      << ss_boundaries[ ii ].second << std::endl;
//        }
//        return ss_boundaries;
//    }
//
//
//    // FT1: DSSP String to FoldTree
//
//    protocols::bootcamp::fold_tree_from_dssp_string(std::string input_dssp_string_) {
//        // 4N - 2 peptide edges, 2N - 2 jump edges
//        utility::vector1< std::pair< core::Size, core::Size > > ss_spans_ = identify_secondary_structure_spans( input_dssp_string_);
//
//        core::Size jedge_count = 1;
//        // first element: jump to mid of every other SS
//        core::kinematics::FoldTree ft;
//        int mid_first_ss = ( ss_spans_[ 1].first + ss_spans_[ 1].second) / 2;
//        for( core::Size i = 2; i <= ss_spans_.size(); ++i) {
//            core::Size mid_i_ss = ( ss_spans_[ i].first + ss_spans_[ i].second) / 2;
//            ft.add_edge( mid_first_ss, mid_i_ss, jedge_count);
//            ++jedge_count;
//        }
//
//        // jumps: mid of first element to mid of every gap
//        for( core::Size i = 1; i < ss_spans_.size(); ++i) {
//            core::Size mid_gap_ss = (ss_spans_[ i].second + ss_spans_[ i+1].first) / 2;
//            ft.add_edge( mid_first_ss, mid_gap_ss, jedge_count);
//            ++jedge_count;
//        }
//
//        for( core::Size i = 1; i <= ss_spans_.size(); ++i) {
//            // each SS: peptide edge mid to end of that SS, peptide edge mid to beginning of that SS
//            core::Size mid_i_ss = ( ss_spans_[ i].first + ss_spans_[ i].second) / 2;
//            ft.add_edge( mid_i_ss, ss_spans_[ i].second, core::kinematics::Edge::PEPTIDE );
//            ft.add_edge( mid_i_ss, ss_spans_[ i].first, core::kinematics::Edge::PEPTIDE );
//        }
//
//        for( core::Size i = 1; i < ss_spans_.size(); ++i) {
//            core::Size mid_gap_ss = (ss_spans_[ i].second + ss_spans_[ i+1].first) / 2;
//            ft.add_edge( mid_gap_ss, ss_spans_[ i].second, core::kinematics::Edge::PEPTIDE);
//            ft.add_edge( mid_gap_ss, ss_spans_[ i+1].first, core::kinematics::Edge::PEPTIDE);
//
//        }
//        ft.add_edge( ss_spans_[ 1].first, 1, core::kinematics::Edge::PEPTIDE);
//        return ft;
//    }
//
//
//    // FT2: SS to String to Fold Tree
//    protocols::bootcamp::fold_tree_from_ss( core::pose::Pose& input_pose_) {
//        core::scoring::dssp::Dssp pose_DSSP( input_pose_);
//        std::string input_dssp_string_ = pose_DSSP.get_dssp_secstruct();
//        return fold_tree_from_dssp_string( input_dssp_string_);
//    }

    void test_fold_tree() {
        std::string input_string = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ";
        core::kinematics::FoldTree ft_test_ = protocols::bootcamp::fold_tree_from_dssp_string( input_string);
        ft_test_.show(std::cout);
        TS_ASSERT_EQUALS( ft_test_.get_jump_edges().size(), 12);
    }

//    void test_fold_tree_Pose() {
//        core::pose::Pose pose_test = create_test_in_pdb_pose();
//        core::kinematics::FoldTree ft_test = fold_tree_from_ss( pose_test);
//        ft_test.show( std::cout);
//        TS_ASSERT( ft_test.check_fold_tree() );
//    }

    void test_hello_world() {
        TS_ASSERT( true);
    }

//    void test_ss_pred_1() {
//        std::string input_string = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
//        utility::vector1< std::pair< core::Size, core::Size > > expected_output =
//                {
//                {4,8},
//                {12,19},
//                {22,26},
//                {36,41},
//                {45,55},
//                {58,62},
//                {65,68}
//                };
//        utility::vector1< std::pair< core::Size, core::Size > > observed_output = identify_secondary_structure_spans( input_string);
//        TS_ASSERT_EQUALS(expected_output, observed_output);
//    }
//
//    void test_ss_pred_2() {
//        std::string input_string = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ";
//        utility::vector1< std::pair< core::Size, core::Size > > expected_output =
//                {
//                        {1,7},
//                        {11,22},
//                        {29,40},
//                        {41,50},
//                        {51,57},
//                        {59,62},
//                        {63,65}
//                };
//        utility::vector1< std::pair< core::Size, core::Size > > observed_output = identify_secondary_structure_spans( input_string);
//        TS_ASSERT_EQUALS(expected_output, observed_output);
//    }
//
//    void test_ss_pred_3() {
//        std::string input_string = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
//        utility::vector1< std::pair< core::Size, core::Size > > expected_output =
//                {
//                        {1,9},
//                        {11,18},
//                        {20,28},
//                        {30,30},
//                        {32,36},
//                        {38,38},
//                        {40,40},
//                        {42,42},
//                        {44,51}
//                };
//        utility::vector1< std::pair< core::Size, core::Size > > observed_output = identify_secondary_structure_spans( input_string);
//        TS_ASSERT_EQUALS(expected_output, observed_output);
//    }

};
