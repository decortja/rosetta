// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyRequirementSet.hh
///
/// @brief A container for all LEGACY_SEWING requirements
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementSet_hh
#define INCLUDED_protocols_legacy_sewing_sampling_requirements_LegacyRequirementSet_hh

//Unit headers
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementSet.fwd.hh>
#include <utility/VirtualBase.hh>
#include <core/types.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.fwd.hh>

#include <protocols/legacy_sewing/sampling/requirements/LegacyGlobalRequirement.fwd.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyIntraSegmentRequirement.fwd.hh>

//Utility headers
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace legacy_sewing  {
namespace sampling {
namespace requirements {

class LegacyRequirementSet : public utility::VirtualBase {

public:

	typedef std::map< core::Size, utility::vector1<LegacyIntraSegmentRequirementOP> > LegacyIntraSegmentRequirementsMap;
	//typedef std::map< std::set<core::Size>, utility::vector1<InterSegmentRequirementOP> > InterSegmentRequirementsMap;

	LegacyRequirementSet();

	core::Size
	min_segments() const;

	void
	min_segments(core::Size min_segments);

	core::Size
	max_segments() const;

	void
	max_segments(core::Size max_segments);

	void
	add_requirement(
		LegacyGlobalRequirementOP requirement
	);

	void
	add_requirement(
		core::Size index,
		LegacyIntraSegmentRequirementOP requirement
	);

	// void
	// add_inter_segment_requirement(
	//  std::set<core::Size> segment
	//  InterSegmentRequirementOP inter_segment_requirement
	// );

	///@brief Evaluated if this Assembly satisfies all
	///contained AssemblyRequirements
	bool satisfies(
		AssemblyCOP assembly
	) const;

	///@brief Evaluate if this Assembly violates any contained
	///AssemblyRequirements
	bool violates(
		AssemblyCOP assembly
	) const;

	///@brief Check all LegacyGlobal requirements to see if we
	///can add more edges to this Assembly
	bool can_be_added_to(
		AssemblyCOP assembly
	) const;

	core::Size get_max_segments()
	const;

	void
	show(
		std::ostream & out
	) const;

private:

	core::Size min_segments_;
	core::Size max_segments_;

	//Requirements for the entire Assembly
	utility::vector1<LegacyGlobalRequirementOP> global_requirements_;

	//Requirements for individual segments
	LegacyIntraSegmentRequirementsMap intra_segment_requirements_;

	//Requirements for interactions between segments
	// utility::vector1< utility::vector <InterSegmentRequirementOP> > inter_segment_requirements_;

};


} //requirements namespace
} //sampling namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
