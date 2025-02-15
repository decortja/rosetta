// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./protocols/pose_creation/MakePolyXMover.hh
/// @brief  header file of MakePolyXMover.cc
/// @author Nobuyasu Koga ( nobuyasau@uw.edu )

#ifndef INCLUDED_protocols_pose_creation_MakePolyXMover_HH
#define INCLUDED_protocols_pose_creation_MakePolyXMover_HH

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/pose_creation/MakePolyXMover.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace pose_creation {


class MakePolyXMover : public protocols::moves::Mover {
public:

	typedef protocols::moves::MoverOP MoverOP;
	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

public:


	// @brief default constructor
	MakePolyXMover();

	// @brief value constructor
	MakePolyXMover( std::string const & aa, bool keep_pro, bool keep_gly, bool keep_disulfide_cys );

	// @brief destructor
	~MakePolyXMover() override;

	/// @brief clone this object
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of object
	protocols::moves::MoverOP fresh_instance() const override;

	// @brief virtual main operation
	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_chis( utility::vector1< core::Real > chis ) { chis_ = chis; }



private:

	/// @brief using amino acid for converting pose to poly XXX
	std::string aa_;

	/// @brief if true, proline, proline are not converted
	bool keep_pro_;

	/// @brief if true, proline, glycine are not converted
	bool keep_gly_;

	/// @brief if true, proline, cystein are not converted
	bool keep_disulfide_cys_;

	/// @brief selects which residues to mutate
	ResidueSelectorCOP selector_;

	/// @brief Use specific chis after mutating
	utility::vector1< core::Real > chis_;

};

} // pose_creation
} // protocols


#endif //INCLUDED_protocols_pose_creation_MakePolyXMover_HH
