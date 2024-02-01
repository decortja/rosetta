// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//

// using the options flags to import a PDB

#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh> 
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>                        // Needed?
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/bootcamp/BootCampMover.hh>


protocols::bootcamp::BootCampMoverOP boot_camp_mover( new protocols::bootcamp::BootCampMover());
protocols::jd2::JobDistributor::get_instance()->go( boot_camp_mover);


//int main ( int argc, char ** argv) {
//
//	// Store a protein's filename and pose.
//	devel::init( argc, argv);
//	utility::vector1< std::string> filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
//	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[ 1]);
//
////    Establishing a bootcamp FoldTree for mypose.
//    mypose->fold_tree( protocols::bootcamp::fold_tree_from_ss( *mypose));
//    mypose->fold_tree().show(std::cout);
//    core::pose::correctly_add_cutpoint_variants( *mypose);
//
//
//	if ( filenames.size() > 0) {
//		std::cout << "You entered: " << filenames[ 1] << " as the PDB file to be read" << std::endl;
//	}
//	else {
//		std::cout << "You didn't provide a PDB file with the -in::file::s option" << std::endl;
//		return 1;
//	}
//
//
//	// Scoring:
//	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
//    sfxn->set_weight( core::scoring::linear_chainbreak, 1);
//	core::Real score = sfxn->score( *mypose);				// mypose is a PoseOP
//	std::cout << "Score: " << score << std::endl;
//
//
//
//    ///////////////////////////////////////////////////////////////////////////
//	// MC: 1 - generate random numbers; 2 - perturb phi/psi; 3 - MC evaluation
//	// MC3: MC object looping
//	protocols::moves::MonteCarlo MC_object_( *mypose, *sfxn, 0.6);
//
//	// MC_ Start pymol
//	protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
//	the_observer->pymol().apply( *mypose);
//
//    // Defining the move map
//	core::kinematics::MoveMap mm;
//	mm.set_bb( true );
//	mm.set_chi( true );
//	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
//	core::optimization::AtomTreeMinimizer atm;
//
//	// Create Copy of pose for speedup (avoids creating/destroying the PoseOP)
//	core::pose::Pose copy_pose = *mypose;
//
//    // Acceptance ratio counter
//    int num_accepted_poses = 0;
//    int num_MC_attempts = 100;
//    core::Real total_score_sum = 0;
//
//	for (int i = 0; i < num_MC_attempts; ++i ) {
//		// MC1: residue selection
//		core::Real uniform_random_number = numeric::random::uniform();
//		core::Size randres = static_cast< core::Size> ( uniform_random_number * mypose->total_residue() + 1);
//
//		// MC2: adjust phi/psi
//		core::Real pert1 = numeric::random::gaussian();
//		core::Real pert2 = numeric::random::gaussian();
//		core::Real orig_phi = mypose->phi( randres);
//		core::Real orig_psi = mypose->psi( randres);
//        mypose->set_phi( randres, orig_phi + pert1);
//		mypose->set_psi( randres, orig_psi + pert2);
//
//        // Repack
//        core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
//        repack_task->restrict_to_repacking();
//        core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
//
//        // Minimize
//        copy_pose = *mypose;
//        atm.run( copy_pose, mm, *sfxn, min_opts );
//        *mypose = copy_pose;
//
//        // Accept or reject, and store counter of accepted poses
//            MC_object_.boltzmann( *mypose);
//            if( MC_object_.boltzmann( *mypose)) {
//                ++num_accepted_poses;
//            }
//            total_score_sum+=MC_object_.last_score();
//	}
//    core::Real fraction_MC_accepted = num_accepted_poses / num_MC_attempts;
//
//    std::cout << "Percent Accepted MC attempts: " << fraction_MC_accepted * 100 << "%" << std::endl;
//	std::cout << "Average score: " << total_score_sum / num_MC_attempts << std::endl;
//
//    return 0;
//}


