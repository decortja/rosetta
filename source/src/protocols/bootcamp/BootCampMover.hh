// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.hh
/// @brief Lab 6: A mover for the monte carlo peptide sampler built in Labs 2 and 4.
/// @author decortja (joseph.a.decorte@vanderbilt.edu)

#ifndef INCLUDED_protocols_bootcamp_BootCampMover_HH
#define INCLUDED_protocols_bootcamp_BootCampMover_HH

// Unit headers
#include <protocols/bootcamp/BootCampMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rosetta_scripts/util.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

// PAREDOWN THIS: TODO
// Unit headers
#include <protocols/bootcamp/BootCampMoverCreator.hh>

#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/PyMOLMover.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// Pare down this list:
#include <devel/init.hh>
#include <numeric/random/random.hh>


namespace protocols {
namespace bootcamp {

///@brief Lab 6: A mover for the monte carlo peptide sampler built in Labs 2 and 4.
class BootCampMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	BootCampMover();
//    BootCampMover( core::Size n_it, core::scoring::ScoreFunctionOP sf );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~BootCampMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap& datamap ) override;

    /// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
    void
    parse_score_function(
            utility::tag::TagCOP const tag,
            basic::datacache::DataMap const & datamap
    );
	//BootCampMover & operator=( BootCampMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

    // Getters and setters for iterations, variables:
    core::Size get_num_iterations() const;
    void set_num_iterations( core::Size num_iterations_new_);

    core::scoring::ScoreFunctionOP get_sfxn() const;
    void set_sfxn( core::scoring::ScoreFunctionOP sfxn_new_);

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns decortja as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

private: // methods

private: // data
    core::scoring::ScoreFunctionOP sfxn_;
    core::Size num_iterations_;
};

std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover );

} //bootcamp
} //protocols

#endif //protocols_bootcamp_BootCampMover_HH
