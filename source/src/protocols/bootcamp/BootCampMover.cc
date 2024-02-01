// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief Lab 6: A mover for the monte carlo peptide sampler built in Labs 2 and 4.
/// @author decortja (joseph.a.decorte@vanderbilt.edu)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>
#include <protocols/jd2/JobDistributor.hh>


// Core headers
#include <core/pose/Pose.hh>
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

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// Pare down this list:
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <utility/pointer/owning_ptr.hh>

// The header for this class
#include <protocols/bootcamp/BootCampMover.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& mypose){

//    Establishing a bootcamp FoldTree for mypose.
            mypose.fold_tree( protocols::bootcamp::fold_tree_from_ss( mypose));
            mypose.fold_tree().show(std::cout);
            core::pose::correctly_add_cutpoint_variants( mypose);

            // Scoring:
            core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
            sfxn->set_weight( core::scoring::linear_chainbreak, 1);
            core::Real score = sfxn->score( mypose);				// mypose is a PoseOP
            std::cout << "Score: " << score << std::endl;

            ///////////////////////////////////////////////////////////////////////////
            // MC: 1 - generate random numbers; 2 - perturb phi/psi; 3 - MC evaluation
            // MC3: MC object looping
            protocols::moves::MonteCarlo MC_object_( mypose, *sfxn, 0.6);

            // MC_ Start pymol
            protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( mypose, true, 0 );
            the_observer->pymol().apply( mypose);

            // Defining the move map
            core::kinematics::MoveMap mm;
            mm.set_bb( true );
            mm.set_chi( true );
            core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
            core::optimization::AtomTreeMinimizer atm;

            // Create Copy of pose for speedup (avoids creating/destroying the PoseOP)
            core::pose::Pose copy_pose = mypose;

            // Acceptance ratio counter
            int num_accepted_poses = 0;
            int num_MC_attempts = 100;
            core::Real total_score_sum = 0;

            for (int i = 0; i < num_MC_attempts; ++i ) {
                // MC1: residue selection
                core::Real uniform_random_number = numeric::random::uniform();
                core::Size randres = static_cast< core::Size> ( uniform_random_number * mypose.total_residue() + 1);

                // MC2: adjust phi/psi
                core::Real pert1 = numeric::random::gaussian();
                core::Real pert2 = numeric::random::gaussian();
                core::Real orig_phi = mypose.phi( randres);
                core::Real orig_psi = mypose.psi( randres);
                mypose.set_phi( randres, orig_phi + pert1);
                mypose.set_psi( randres, orig_psi + pert2);

                // Repack
                core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( mypose );
                repack_task->restrict_to_repacking();
                core::pack::pack_rotamers( mypose, *sfxn, repack_task );

                // Minimize
                copy_pose = mypose;
                atm.run( copy_pose, mm, *sfxn, min_opts );
                mypose = copy_pose;

                // Accept or reject, and store counter of accepted poses
                MC_object_.boltzmann( mypose);
                if( MC_object_.boltzmann( mypose)) {
                    ++num_accepted_poses;
                }
                total_score_sum+=MC_object_.last_score();
            }
            core::Real fraction_MC_accepted = num_accepted_poses / num_MC_attempts;

            std::cout << "Percent Accepted MC attempts: " << fraction_MC_accepted * 100 << "%" << std::endl;
            std::cout << "Average score: " << total_score_sum / num_MC_attempts << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
    /* attlist + XMLSchemaAttribute( "niterations", core::Size, "The number of Monte Carlo iterations to run.")
    + attribute_for_parse_score_function; */ ////TODO
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Lab 6: A mover for the monte carlo peptide sampler built in Labs 2 and 4.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );

}

/// @brief This mover is unpublished.  It returns decortja as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"decortja",
		"TODO: institution",
		"joseph.a.decorte@vanderbilt.edu",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
