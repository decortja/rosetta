// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hbnet/hbnet_util.hh
/// @brief utility functions used by many of the HBNet movers and filters
/// @author Scott Boyken (sboyken@gmail.com)


#ifndef INCLUDED_protocols_hbnet_HBNet_util_hh
#define INCLUDED_protocols_hbnet_HBNet_util_hh

#include <protocols/hbnet/HBNet.fwd.hh>

#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>

#include <core/types.hh>
//#include <core/graph/Graph.hh>
#include <core/scoring/hbonds/graph/HBondGraph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>

#include <utility/cxx_versioning_macros.hh>

#include <protocols/hbnet/NetworkState.fwd.hh> // AUTO IWYU For NetworkState
#include <core/pack/rotamer_set/RotamerSet.fwd.hh> // AUTO IWYU For RotamerSet
#include <core/pack/task/PackerTask.fwd.hh> // AUTO IWYU For PackerTask
#include <core/scoring/ScoreFunction.hh> // AUTO IWYU For ScoreFunction

namespace protocols {
namespace hbnet {

PURE_ATTR
std::string
print_list_to_string(
	HBondNetStruct const & network,
	bool chainid=true,
	bool term_w_start=false,
	bool term_w_cycle=false,
	bool term_w_bb=false
);

PURE_ATTR
std::string print_list_to_string(
	core::pose::Pose const & pose,
	HBondNetStruct const & network,
	bool chainid=true,
	bool term_w_start=false,
	bool term_w_cycle=false,
	bool term_w_bb=false,
	bool use_pdb_numbering=true
);

PURE_ATTR
std::string print_network( HBondNetStruct const & i, bool chainid=true );

PURE_ATTR
std::string print_network_w_pdb_numbering( core::pose::Pose const & pose, HBondNetStruct const & i, bool chainid );

CONST_ATTR NODISCARD_ATTR
std::string print_headers();

PURE_ATTR
utility::vector1< HBondResStructCOP >::const_iterator
find_HBondResStruct(
	utility::vector1< HBondResStructCOP > const & residues,
	core::Size resnum
);

PURE_ATTR
void get_hbond_atom_pairs(
	HBondNetStruct & network,
	core::pose::Pose & pose,
	bool bb_exlcusion=false
);

PURE_ATTR
bool hbond_exists_in_vector(
	utility::vector1<core::scoring::hbonds::HBondCOP> const & hbond_vec,
	core::scoring::hbonds::HBondCOP & h2
);

PURE_ATTR
void add_reslabels_to_pose(
	core::pose::Pose & pose,
	HBondNetStruct & i,
	std::string label="HBNet"
);

PURE_ATTR
core::Size get_num_protein_sc_sc_hbonds(
	core::pose::Pose & pose,
	HBondNetStruct & i
);

PURE_ATTR
core::Size get_num_edges_for_res(
	core::Size const res,
	ObjexxFCL::FArray2D_int & path_dists
);

PURE_ATTR
void hbnet_symm_one_body_energies(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSet & rotset_op,
	core::scoring::ScoreFunction const & sf,
	core::pack::task::PackerTask const & task,
	utility::graph::Graph const & packer_neighbor_graph,
	utility::vector1< core::PackerEnergy > & energies
);

PURE_ATTR
void hbnet_one_body_energies(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSet & rotset_op,
	core::scoring::ScoreFunction const & sf,
	utility::vector1< core::PackerEnergy > & energies
);

PURE_ATTR
bool network_contains_aa(
	char aa_one_letter,
	HBondNetStruct const & i
);

PURE_ATTR
bool network_contains_aa(
	char aa_one_letter,
	utility::vector1< HBondResStructCOP > const & residues
);

PURE_ATTR
bool his_tyr_connectivity( core::pose::Pose const & pose, HBondNetStruct & i );


PURE_ATTR
bool
edge_satisfies_heavy_unsat_for_node(
	NetworkState const & current_state,
	core::scoring::hbonds::graph::HBondNode const * node,
	core::scoring::hbonds::graph::HBondEdge const * edge
);

} //hbnet
} //protocols


#endif
