// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/MergePDBMover.hh
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_protocols_pose_creation_MergePDBMover_hh
#define INCLUDED_protocols_pose_creation_MergePDBMover_hh

#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/pose_creation/MergePDBMover.fwd.hh>
#include <utility/vector1.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace pose_creation {

class MergePDBMover : public moves::Mover {

public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSelectorCOP ResidueSelectorCOP;

	struct Overlap{
		core::Size start_overlap_xmlPose;
		core::Size end_overlap_xmlPose;
		core::Size start_overlap_pose;
		core::Size end_overlap_pose;
		core::pose::PoseOP output_poseOP;
		bool output_yet;
		Overlap(core::Size start_overlap_xmlPose_i, core::Size end_overlap_xmlPose_i, core::Size start_overlap_pose_i, core::Size end_overlap_pose_i){
			start_overlap_xmlPose = start_overlap_xmlPose_i;
			end_overlap_xmlPose = end_overlap_xmlPose_i;
			start_overlap_pose = start_overlap_pose_i;
			end_overlap_pose = end_overlap_pose_i;
			output_yet=false;
			output_poseOP=nullptr;
		}

		std::string name(){

			std::string name_str = std::to_string(start_overlap_pose) + "_" + std::to_string(end_overlap_pose) + "_" + std::to_string(start_overlap_xmlPose) + "_" + std::to_string(end_overlap_xmlPose);
			return(name_str);
		}

		void print(){
			std::cout << start_overlap_xmlPose <<"," << end_overlap_xmlPose << "::" << start_overlap_pose << "," << start_overlap_pose << std::endl;
		}
	};
	/// @brief  Constructor
	MergePDBMover();
	/// @brief Determines the overlaps. stores the start and end position in the struct Overlap
	void determine_overlap(Pose const pose,core::Size chain_id);
	/// @brief fast check to maake sure the pose isn't a structural duplicate
	bool check_duplicate(Pose & pose);
	/// @brief Uses the overlap to generate poses
	void generate_overlaps(Pose & pose,core::Size chain_id);
	/// @brief Figures out the closest residue that's not part of the overlap.
	core::Size closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid);
	/// @brief Gets the entire SS element the match is on
	void increase_range_to_ignore_ss_element(core::pose::Pose const & pose, core::Size init_start, core::Size init_end, core::Size & ss_start, core::Size & ss_end);
	/// @brief Copies the sequence in the overlap region as appropriate. This sets the initial residues before the pack and minimize is called
	void copy_sequence(core::Size start_overlap_resid, core::Size end_overlap_resid, core::Size start_overlap_input_pose_resid,core::Size start_overlap_xml_pose_resid,core::pose::Pose const & input_pose,core::pose::Pose const & xml_pose,core::pose::Pose & output_pose);
	/// @brief packs and minimizes if no clashes as determined by score0
	void pack_and_minimize(Pose const pose, core::Real baseline_score);
	/// @brief any poses that score lower than the input files is output
	core::pose::PoseOP get_additional_output() override;

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	utility::vector1<MergePDBMover::Overlap> overlaps_;
	core::pose::PoseOP xml_input_pose_;
	std::string overlap_location_pose_;
	core::Real overlap_max_rmsd_;
	core::Size overlap_length_;
	core::Size overlap_scan_range_cmdLine_;
	core::Size overlap_scan_range_xml_;
	core::Real design_range_;
	core::Real packing_range_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::scoring::ScoreFunctionOP asymm_score_;
	core::pack::task::TaskFactoryOP task_factory_;
	bool do_minimize_;
	std::string chain_;
	std::string symm_file_;
	std::string no_design_label_;
	std::string init_overlap_sequence_;
	core::Real duplicate_rmsd_pose_threshold_;
	bool detect_disulf_before_repack_;
	bool output_only_first_;
	bool output_overlap_positions_;
	bool do_design_;
	core::Real clash_threshold_;
	core::select::residue_selector::ResidueSelectorCOP selector_cmd_;
	core::select::residue_selector::ResidueSelectorCOP selector_xml_;
};

} // pose_creation
} // protocols

#endif

