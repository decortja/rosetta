// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/bootcamp/fold_tree_from_ss.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

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
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// C++ headers
#include <string>


namespace protocols {
    namespace bootcamp {
        // FT1: DSSP String to FoldTree
        utility::vector1< std::pair< core::Size, core::Size > > identify_secondary_structure_spans( std::string const & ss_string )
        {
            utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
            core::Size strand_start = -1;
            for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
                if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
                    if ( int( strand_start ) == -1 ) {
                        strand_start = ii;
                    } else if ( ss_string[ii] != ss_string[strand_start] ) {
                        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                        strand_start = ii;
                    }
                } else {
                    if ( int( strand_start ) != -1 ) {
                        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                        strand_start = -1;
                    }
                }
            }
            if ( int( strand_start ) != -1 ) {
                // last residue was part of a ss-element
                ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
            }
            for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
                std::cout << "SS Element " << ii << " from residue "
                          << ss_boundaries[ ii ].first << " to "
                          << ss_boundaries[ ii ].second << std::endl;
            }
            return ss_boundaries;
        }
        core::kinematics::FoldTree fold_tree_from_dssp_string(std::string input_dssp_string_) {
            // 4N - 2 peptide edges, 2N - 2 jump edges
            utility::vector1< std::pair< core::Size, core::Size > > ss_spans_ = identify_secondary_structure_spans( input_dssp_string_);

            core::Size jedge_count = 1;
            core::kinematics::FoldTree ft;
            core::Size mid_first_ss = floor( ss_spans_[ 1].first + ss_spans_[ 1].second) / 2;


            // Jump Edges:
            // jump edge: mid first element to mid of every other SS
            for( core::Size i = 2; i <= ss_spans_.size(); ++i) {
                core::Size mid_i_ss = floor( (ss_spans_[ i].first + ss_spans_[ i].second) / 2);
                ft.add_edge( mid_first_ss, mid_i_ss, jedge_count);
                ++jedge_count;
            }

            // jump edge: mid of first element to mid of every gap
            for( core::Size i = 1; i < ss_spans_.size(); ++i) {
                // make sure a gap actually exists
                if( ( ss_spans_[ i+1].first - ss_spans_[ i].second) > 1) {
                    core::Size mid_gap_ss = (ss_spans_[ i].second + ss_spans_[ i+1].first) / 2;
                    ft.add_edge( mid_first_ss, mid_gap_ss, jedge_count);
                    ++jedge_count;
                }
                else {
                    continue;
                }
            }

            // Peptide Edges:
            // peptide edge: mid of first element to first residue
            core::Size begin = 1;

            // peptide edge: mid of every SS to end of that SS, mid of every SS to beginning of that SS
            for( core::Size i = 1; i <= ss_spans_.size(); ++i) {
                // Left-direction:
                core::Size mid_i_ss = (ss_spans_[i].first + ss_spans_[i].second) / 2;
                if( i == 1) {
                    ft.add_edge( mid_i_ss, begin, core::kinematics::Edge::PEPTIDE);
                }
                else {
                    ft.add_edge(mid_i_ss, ss_spans_[i].first, core::kinematics::Edge::PEPTIDE);
                }

                // Right direction:
                if( i == ss_spans_.size()) {
                    ft.add_edge( mid_i_ss, input_dssp_string_.length(), core::kinematics::Edge::PEPTIDE);
                }
                else {
                    ft.add_edge(mid_i_ss, ss_spans_[i].second, core::kinematics::Edge::PEPTIDE);
                }
            }

            // For gaps:
            for( core::Size i = 1; i < ss_spans_.size(); ++i) {
                if( (ss_spans_[ i+1].first - ss_spans_[ i].second) > 1) {
                    core::Size mid_gap_ss = (ss_spans_[ i].second + ss_spans_[ i+1].first) / 2;
                    ft.add_edge( mid_gap_ss, ss_spans_[ i].second + 1, core::kinematics::Edge::PEPTIDE);
                    ft.add_edge( mid_gap_ss, ss_spans_[ i+1].first - 1, core::kinematics::Edge::PEPTIDE);
                }
            }
            return ft;
        }


        // FT2: SS to String to Fold Tree
        core::kinematics::FoldTree fold_tree_from_ss( core::pose::Pose& input_pose_) {
            core::scoring::dssp::Dssp pose_DSSP( input_pose_);
            std::string input_dssp_string_ = pose_DSSP.get_dssp_secstruct();
            return fold_tree_from_dssp_string( input_dssp_string_);
        }
    }
}