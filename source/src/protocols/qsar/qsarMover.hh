// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/qsar/qsarMover.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_qsarMover_hh
#define INCLUDED_protocols_qsar_qsarMover_hh

#include <protocols/qsar/qsarMover.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>


#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace qsar {

class qsarMover : public protocols::moves::Mover
{
public:
	qsarMover();
	//qsarMover(scoring_grid::GridSetCOP grid_set_prototype, core::Real width, core::Real resolution);
	void apply(core::pose::Pose & pose) override;
	void parse_my_tag(utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;
	std::string get_name() const override;

private:
	scoring_grid::GridSetCOP grid_set_prototype_;
	qsarMapOP qsar_map_;
	std::string chain_;
	utility::vector1<std::string>  grids_to_use_;
};

}
}

#endif /* QSARMOVER_HH_ */
