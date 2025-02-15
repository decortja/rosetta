// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
/// @brief demo program for implementing loop relax + FA relax
/// @author Mike Tyka
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Daniel J. Mandell

// include these first for building on Visual Studio

#include <protocols/relax/loop/LoopRelaxMover.hh>
#include <protocols/relax/loop/LoopRelaxMoverCreator.hh>

// Package headers
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/chemical/VariantType.hh>
#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif // BOINC_GRAPHICS

// Project Headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/tag/Tag.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/FragSet.hh>
#include <core/id/AtomID.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/loops/loops_definers/util.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/LoopMoverFactory.hh>
//#include <protocols/loops/LoopMover.fwd.hh>
//#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD.hh>
#include <protocols/loops/looprelax_protocols.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// Refactored Kic headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/samplers/KicSampler.hh>
#include <protocols/loop_modeling/samplers/BalancedKicSampler.hh>
#include <protocols/loop_modeling/samplers/LegacyKicSampler.hh>
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/utilities/RepeatedTask.hh>
#include <protocols/loop_modeling/utilities/PeriodicTask.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/ProgressBar.hh>
#include <protocols/loop_modeling/loggers/AcceptanceRates.hh>
#include <protocols/loop_modeling/loggers/PdbLogger.hh>
#include <protocols/loop_modeling/loggers/ScoreVsRmsd.hh>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>  // this was auto-removed but is needed for graphics builds!
#endif
#include <utility/exit.hh>

// symmetry
#include <core/pose/symmetry/util.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>


// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/init_id_map.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/util.hh>
#include <utility/Bound.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/options/BooleanOption.hh>
#include <basic/options/option.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif
// }}}1

namespace protocols {
namespace comparative_modeling {

// {{{1
LoopRelaxMover::LoopRelaxMover() : moves::Mover(),
	guarded_loops_( new loops::GuardedLoopsFromFile )
{
	set_defaults_();
}

// {{{1
/// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
/// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
/// NOT BE INITIALIZED
LoopRelaxMover::LoopRelaxMover(
	std::string const & remodel,
	std::string const & intermedrelax,
	std::string const & refine,
	std::string const & relax,
	loops::Loops const & loops
) : moves::Mover(),
	cmd_line_csts_( true ),
	copy_sidechains_( true ),
	n_rebuild_tries_( 3 ),
	rebuild_filter_( 999 ),
	remodel_( remodel ),
	intermedrelax_( intermedrelax ),
	refine_( refine ),
	relax_( relax ),
	guarded_loops_( new loops::GuardedLoopsFromFile( loops ))
{}

// {{{1
/// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
/// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
/// NOT BE INITIALIZED
LoopRelaxMover::LoopRelaxMover(
	std::string const & remodel,
	std::string const & intermedrelax,
	std::string const & refine,
	std::string const & relax,
	loops::LoopsFileData const & loops_from_file
) : moves::Mover(),
	cmd_line_csts_( true ),
	copy_sidechains_( true ),
	n_rebuild_tries_( 3 ),
	rebuild_filter_( 999 ),
	remodel_( remodel ),
	intermedrelax_( intermedrelax ),
	refine_( refine ),
	relax_( relax ),
	guarded_loops_( new loops::GuardedLoopsFromFile( loops_from_file ))
{}

// {{{1
LoopRelaxMover::LoopRelaxMover(
	std::string const & remodel,
	std::string const & intermedrelax,
	std::string const & refine,
	std::string const & relax,
	loops::GuardedLoopsFromFileOP guarded_loops
) : moves::Mover(),
	cmd_line_csts_( true ),
	copy_sidechains_( true ),
	n_rebuild_tries_( 3 ),
	rebuild_filter_( 999 ),
	remodel_( remodel ),
	intermedrelax_( intermedrelax ),
	refine_( refine ),
	relax_( relax ),
	guarded_loops_( guarded_loops ) // shallow copy
{}

// {{{1
/// @brief Copy-ctor; shallow copy of all data object.
LoopRelaxMover::LoopRelaxMover( LoopRelaxMover const & src ) : moves::Mover(),
	cmd_line_csts_( src.cmd_line_csts_ ),
	copy_sidechains_( src.copy_sidechains_ ),
	n_rebuild_tries_( src.n_rebuild_tries_ ),
	rebuild_filter_( src.rebuild_filter_ ),
	remodel_( src.remodel_ ),
	intermedrelax_( src.intermedrelax_ ),
	refine_( src.refine_ ),
	relax_( src.relax_ ),
	guarded_loops_( src.guarded_loops_ ),
	cen_scorefxn_( src.cen_scorefxn_ ),
	fa_scorefxn_( src.fa_scorefxn_ ),
	frag_libs_( src.frag_libs_ ),
	compute_rmsd_( src.compute_rmsd_ )
{}

// {{{1
/// @brief assignment operator; Shallow copy of all data.
LoopRelaxMover const & LoopRelaxMover::operator=( LoopRelaxMover const & rhs )
{
	if ( this != &rhs ) {
		cmd_line_csts_ = rhs.cmd_line_csts_;
		copy_sidechains_ = rhs.copy_sidechains_;
		n_rebuild_tries_ = rhs.n_rebuild_tries_;
		rebuild_filter_ = rhs.rebuild_filter_;
		remodel_ = rhs.remodel_;
		intermedrelax_ = rhs.intermedrelax_;
		refine_ = rhs.refine_;
		relax_ = rhs.relax_;
		guarded_loops_ = rhs.guarded_loops_;
		cen_scorefxn_ = rhs.cen_scorefxn_;
		fa_scorefxn_ = rhs.fa_scorefxn_;
		frag_libs_ = rhs.frag_libs_;
		compute_rmsd_ = rhs.compute_rmsd_;
	}
	return *this;
}

// {{{1
LoopRelaxMover::~LoopRelaxMover() {}

// {{{1
void LoopRelaxMover::frag_libs(
	utility::vector1< core::fragment::FragSetOP > new_libs
) {
	frag_libs_ = new_libs;
}

// {{{1
utility::vector1< core::fragment::FragSetOP > LoopRelaxMover::frag_libs() const {
	return frag_libs_;
}

// {{{1
void LoopRelaxMover::scorefxns(
	core::scoring::ScoreFunctionOP centroid_scorefxn,
	core::scoring::ScoreFunctionOP fullatom_scorefxn
) {
	//cen_scorefxn_ = cen_scorefxn;
	//fa_scorefxn_ = fa_scorefxn;
	cen_scorefxn( centroid_scorefxn );
	fa_scorefxn ( fullatom_scorefxn );
}

// {{{1
void LoopRelaxMover::fa_scorefxn(
	core::scoring::ScoreFunctionOP fa_scorefxn
) {
	fa_scorefxn_ = fa_scorefxn;
}

// {{{1
void LoopRelaxMover::cen_scorefxn(
	core::scoring::ScoreFunctionOP cen_scorefxn
) {
	cen_scorefxn_ = cen_scorefxn;
}

// {{{1
void LoopRelaxMover::apply( core::pose::Pose & pose ) {
	int corelength = 0;
	std::string const curr_job_tag( get_current_tag() );

	// store the initial (possibly fullatom) pose
	//   used for computing RMS later
	//   we may steal sidechains from this pose as well
	core::pose::Pose start_pose = pose;

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// this typedef belongs in its own .fwd.hh file.
	//typedef utility::pointer::owning_ptr< loops::IndependentLoopMover > loops::IndependentLoopMoverOP;

	basic::Tracer TR( "protocols.looprelax" );

	TR << "==== Loop protocol: ================================================="
		<< std::endl;
	TR << " remodel        " << remodel()        << std::endl;
	TR << " intermedrelax  " << intermedrelax()  << std::endl;
	TR << " refine         " << refine()         << std::endl;
	TR << " relax          " << relax()          << std::endl;

	// load native pose (if provided)
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		core::pose::set_ss_from_phipsi( native_pose );
	} else {
		native_pose = start_pose;
	}

	if ( start_pose.size() != native_pose.size() ) {
		// strip VRTs from the end, then compare lengths
		int nnonvrt_start = start_pose.size();
		while ( nnonvrt_start>0 && start_pose.residue( nnonvrt_start ).aa() == core::chemical::aa_vrt ) nnonvrt_start--;

		int nnonvrt_native = native_pose.size();
		while (  nnonvrt_native>0 && native_pose.residue( nnonvrt_native ).aa() == core::chemical::aa_vrt ) nnonvrt_native--;
		if ( nnonvrt_native != nnonvrt_start ) {
			utility_exit_with_message(
				"Start pose and native pose don't match in length"
			);
		}
	}

	evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
	evaluator->add_evaluation(
		new simple_filters::SelectRmsdEvaluator( native_pose, "_native" )
	);

#ifdef BOINC_GRAPHICS
	// set native for graphics
	boinc::Boinc::set_graphics_native_pose( native_pose );
#endif

	// pick loops if necessary
	guarded_loops_->resolve_loop_indices( pose );
	protocols::loops::LoopsOP loops = get_loops();
	// try to load loops from command line
	if ( loops->size() == 0 ) {
		TR.Debug << "picking loops by chainbreak score." << std::endl;
		protocols::loops::LoopsOP loops_picked = protocols::loops::pick_loops_chainbreak(
			start_pose, option[ cm::min_loop_size ]() );
		*loops = *loops_picked; // copy the loops we just read in into the guarded_loops_ object

		if ( loops->size() == 0 ) {
			TR.Debug << "no loops found." << std::endl;
			remodel( "no" );
		}
	} // loops.size() == 0

	loops->verify_against( start_pose );
	TR.Debug << loops << std::endl;

	if ( option[ OptionKeys::loops::extended ]() ) loops->set_extended( true );
	bool const debug( option[ OptionKeys::loops::debug ]() );

	// superimpose native over core ?
	core::pose::Pose native_pose_super = native_pose;
	id::AtomID_Map< id::AtomID > atom_map;
	if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
		core::pose::initialize_atomid_map( atom_map, native_pose_super, core::id::AtomID::BOGUS_ATOM_ID() );
		for ( core::Size ir=1; ir <= native_pose.size(); ++ir ) {
			if ( !loops->is_loop_residue( ir ) ) {
				id::AtomID const id1( native_pose_super.residue(ir).atom_index("CA"), ir );
				id::AtomID const id2( pose.residue(ir).atom_index("CA"), ir );
				atom_map.set(id1, id2);
			}
		}
		/*core::Real rms =*/ core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_super_pose.pdb");
		if ( debug ) native_pose_super.dump_pdb(curr_job_tag + "_after_super_native.pdb");
	}

	bool const fullatom_input( pose.is_fullatom() );
	bool const fullatom_output(
		option[ out::file::fullatom ]() || refine() != "no" || relax() != "no"
	);

	if ( fullatom_input ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
	}


	checkpoint::CheckPointer checkpoints_("Loopbuild");
	// need to make sure that the bailout structure is set!
	checkpoints_.checkpoint( pose, curr_job_tag, "initial", true );

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

#ifdef GL_GRAPHICS
	protocols::viewer::add_conformation_viewer(
		pose.conformation(), "loops_pose"
	);
#endif
	////////////////////////////

	// loop rebuilding
	loops->auto_choose_cutpoints( pose );

	// if remove_extended_loops is specified, treat extended loops as missing density, by randomly placing the atoms
	// this is not great behavior, BUT the code is already tolerant of missing density treated in this fashion
	bool remove_extended_loops = option[ OptionKeys::loops::remove_extended_loops ]();
	if ( remove_extended_loops ) {
		for ( loops::Loops::const_iterator it = loops->begin(), it_end = loops->end();
				it != it_end; ++it
				) {
			if ( it->is_extended() && it->skip_rate() == 0.0 ) {
				TR << "Removing loop: " << *it << std::endl;
				int lstart = it->start(); if ( lstart != 1 ) lstart++;

				for ( core::Size r = lstart; r<= it->stop(); ++r ) {
					for ( core::Size k = 1; k<= pose.residue(r).natoms(); ++k ) {
						numeric::xyzVector< core::Real > rnd_atm(
							900.000 + numeric::random::uniform()*100.000,
							900.000 + numeric::random::uniform()*100.000,
							900.000 + numeric::random::uniform()*100.000
						);
						pose.set_xyz( core::id::AtomID( k,r ) , rnd_atm );
					}
				}
			}
		}
	} // remove_extended_loops

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////   Remove Missing density
	////
	//// if this is an initial loop building exercise, build conservatively first, i.e. dont extend loop
	//// regions etc. Then repeat with a more aggressive approach to make sure the loops are closed etc.
	//// This makes sure all the missing density is gone before doing proper loop building.
	////
	if ( option[ OptionKeys::loops::build_initial ].user() ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Initial Building     " << std::endl;
		TR << "===" << std::endl;

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "initial_build", false, true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_initial_build.pdb");

			runtime_assert( frag_libs().size() > 0 );
			loops::loop_mover::perturb::LoopMover_Perturb_QuickCCD quick_ccd( loops );
			for ( Size i = 1; i <= frag_libs().size(); ++i ) {
				quick_ccd.add_fragments( frag_libs()[i] );
			}

			quick_ccd.get_checkpoints()->set_type("InitialBuild");
			quick_ccd.set_current_tag( curr_job_tag );
			quick_ccd.set_native_pose( new core::pose::Pose ( native_pose ) );
			quick_ccd.set_scorefxn( cen_scorefxn_ );
			quick_ccd.set_build_attempts_( 1 );
			quick_ccd.set_grow_attempts_( 0 );
			quick_ccd.set_accept_aborted_loops_( true );
			quick_ccd.set_strict_loops( true );
			quick_ccd.set_random_order_( false );
			quick_ccd.set_build_all_loops_( true );
			quick_ccd.set_loop_combine_rate_( 0.0 );

			/// RUN quick_ccd
			quick_ccd.apply( pose );

			loops::remove_cutpoint_variants( pose );

			//fpd if we care at all about the global coordinate frame
			//fpd (and thus have a root VRT) then don't recenter the pose
			if ( pose.residue_type( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
				pose.center();
			}
			(*cen_scorefxn_)(pose);
			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_initial_build.pdb");

			checkpoints_.checkpoint( pose, curr_job_tag, "initial_build", true);
		}
		checkpoints_.debug( curr_job_tag, "initial_build", (*cen_scorefxn_)( pose ) );
	} // build_initial

	// Make sure loops can be grown in any protocol, not just QuickCCD
	if ( basic::options::option[ basic::options::OptionKeys::loops::random_grow_loops_by ].user() ) {
		loops->grow_all_loops( pose ,  basic::options::option[ basic::options::OptionKeys::loops::random_grow_loops_by ]() );
	}

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  add constraints if specified by user.
	////
	////
	if ( cmd_line_csts() ) {
		core::scoring::constraints::add_constraints_from_cmdline(
			pose, *cen_scorefxn_
		);
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(
			*cen_scorefxn_
		);
	}

	// mjo TODO: check if pose.conformation().detect_disulfides needs to be called.  If it does consider replacing this with core/pose/initialize_disulfide_bonds().

	// read in disulfides if specified by user
	if ( option[ in::fix_disulf ].user() ) {
		using std::pair;
		using utility::vector1;
		vector1< pair<Size, Size> > disulfides;

		core::io::raw_data::DisulfideFile ds_file( option[ in::fix_disulf ]() );
		ds_file.disulfides( disulfides, pose );
		pose.conformation().fix_disulfides( disulfides );
	}


	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Loop remodelling (centroid loop modelling)
	////
	////


	long starttime = time(NULL);
	bool all_loops_closed = true;
	bool tmp_all_loops_closed = false;
	if ( remodel() != "no" ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Remodel    " << std::endl;
		TR << "===" << std::endl;
		all_loops_closed = false;
		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "remodel", false, true) ) {

			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_rebuild.pdb");

			using core::Size;
			using core::Real;
			using core::pose::Pose;
			for ( Size ii = 1; ii <= n_rebuild_tries(); ++ii ) {
				core::Real current_sc( rebuild_filter() + 1 );
				TR.Debug << "Remodeling attempt " << ii << "." << std::endl;
				if ( remodel() == "old_loop_relax" ) {
					LoopRebuild loop_rebuild( cen_scorefxn_, *loops );
					loop_rebuild.apply( pose );
				} else {

					/* // DJM: does this cause a crash if the only loop is terminal
					if ( remodel() == "perturb_kic" ) {
					// remove the terminal loops - perturb_kic doesnt seem to support these, sadly.
					protocols::loops::Loops newloops;
					for( core::Size i=1; i <= loops.size() ; i ++ ){
					if( loops[i].start() <= 1 ) continue;
					if( loops[i].stop() >= pose.size() ) continue;
					newloops.add_loop( loops[i] );
					}
					loops = newloops;
					}
					*/
					// DJM: need to cast this as IndependentLoopMover to set strict loops to true.
					loops::loop_mover::IndependentLoopMoverOP remodel_mover( static_cast< loops::loop_mover::IndependentLoopMover * >
						( loops::LoopMoverFactory::get_instance()->create_loop_mover( remodel(), loops ).get() ) );
					core::kinematics::FoldTree f_orig=pose.fold_tree();
					if ( !remodel_mover ) {
						utility_exit_with_message( "Error: no remodel mover defined!" );
					}
					if ( ! ( remodel() == "perturb_kic" ) ) {
						runtime_assert( frag_libs().size() > 0 );
					}

					for ( Size i = 1; i <= frag_libs().size(); ++i ) {
						remodel_mover->add_fragments( frag_libs()[i] );
					}

					remodel_mover->set_scorefxn( cen_scorefxn_ );  // set cst score

					if ( remodel() == "perturb_kic" ) {
						core::kinematics::FoldTree f_new;
						protocols::loops::fold_tree_from_loops( pose, *loops,  f_new, true );
						pose.fold_tree( f_new );
					}
					remodel_mover->get_checkpoints()->set_type("Remodel");
					remodel_mover->set_current_tag( curr_job_tag );
					remodel_mover->set_native_pose( new Pose( native_pose ) );
					remodel_mover->apply( pose );

					if ( remodel() == "perturb_kic" ) { //DJM: skip this struct if initial closure fails
						if ( remodel_mover->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
							set_last_move_status(protocols::moves::FAIL_RETRY);
							TR << "Structure " << " failed initial kinematic closure. Skipping..." << std::endl;
							//bool fail = true; // make this
							pose.fold_tree( f_orig );
							return;
						}
						pose.fold_tree( f_orig );
					}//if ( remodel() == "perturb_kic" )
					tmp_all_loops_closed = remodel_mover->get_all_loops_closed();
				}
				current_sc = (*cen_scorefxn_)( pose );

				if ( option[ OptionKeys::cm::loop_rebuild_filter ].user() || (remodel() == "old_loop_relax") ) {
					TR << "classic check for loop closure" << std::endl;
					current_sc = (*cen_scorefxn_)( pose );
					if ( current_sc <= rebuild_filter() ) {
						all_loops_closed = true;
						break;  //fpd changed >= to <=
						// ... don't we want to stop when our score is _less_ than the cutoff
					}
				} else {
					if ( tmp_all_loops_closed ) {
						all_loops_closed = true;
						break;
					}
				}
			} // for ( in n_rebuild_tries )

			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_rebuild.pdb");
			checkpoints_.checkpoint( pose, curr_job_tag, "remodel", true);
		} // recover checkpoint
		checkpoints_.debug( curr_job_tag, "remodel", (*cen_scorefxn_)( pose ) );
	} else { all_loops_closed = true; } // if ( remodel != no ) // AS 03/16/2012: hack to allow refinement without a preceding perturb stage

	long endtime = time(NULL);

	TR << "Buildtime: " << endtime - starttime << std::endl;


	TR << pose.fold_tree() << std::endl;

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Halfway stats
	////
	////
	if ( compute_rmsd() ) {
		setPoseExtraScore( pose, "cen_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
		if ( option[ in::file::native ].user() ) {
			setPoseExtraScore( pose, "cen_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
			setPoseExtraScore( pose, "cen_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
			setPoseExtraScore( pose, "cen_loopcarms",  loops::loop_rmsd(native_pose_super, pose, *loops, true ) );
		}
	}


	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Full atom part
	////
	////
	if ( fullatom_output ) {
		// make sure there aren't any cut point variants flying about

		loops::remove_cutpoint_variants( pose, true );

		// if no centroid modelling was done at all, then grab the original
		// fullatom pose.
		if ( remodel() == "no" &&
				!option[ OptionKeys::loops::build_initial ].user() &&
				fullatom_input
				) {
			pose = start_pose;
		}

		TR << "===================================================================================="
			<< std::endl;
		TR << "===" << std::endl;
		TR << "===   Fullatom " << std::endl;
		TR << "===" << std::endl;

		if ( debug ) pose.dump_pdb(curr_job_tag + "_before_fullatom.pdb");
		TR << "Annotated sequence before fa switch: " << pose.annotated_sequence(true) << std::endl;

		//puts full-atom sidechains on loop regions
		core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
		pose.conformation().detect_bonds(); //apl fix this !

		utility::vector1< bool > needToRepack( pose.size() , !fullatom_input );
		bool needToRepackAtAll = !fullatom_input;

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_fullatom.pdb");
		// copy sidechain torsions from input pose
		if ( fullatom_input && copy_sidechains() )  {

			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				// if remodelling was done, repack the loops - otherwise leave it.
				if ( remodel() != "no" ) {
					for ( loops::Loops::const_iterator it=loops()->begin(), it_end=loops()->end(); it != it_end; ++it ) {
						if (    i >= core::Size( it->start() ) - 3
								&& i <= core::Size( it->stop() ) + 3 ) {
							// allow 3-residue leeway on either side for 'random_loops'
							// this kind of sucks.
							TR.Debug << "Repacking because in loop: " << i << std::endl;
							needToRepack[i] = true;
							break;
						}
					}
				}

				// if there is missing density in the sidechain, then we need to
				// repack check this my making sure that no SC atom is more than
				// 20A (?) away from CA
				if ( start_pose.residue_type(i).is_protein() && start_pose.residue_type(i).has("CA") ) {
					numeric::xyzVector< core::Real> ca_pos = start_pose.residue(i).atom("CA").xyz();
					for ( int j=1; j<=(int)start_pose.residue(i).natoms(); ++j ) {
						if ( (ca_pos - start_pose.residue(i).atom(j).xyz()).length() > 20 ) {
							TR.Debug << "Missing dens: " << i << std::endl;
							needToRepack[i] = true;
							break;
						}
					}
				} // missing density?

				//copy sidechains only for non-loop regions
				if ( !needToRepack[i] ) {
					if ( pose.residue_type(i).is_protein() ) {
						pose.replace_residue( i, start_pose.residue(i), true );
					}
					TR.Debug << "Copying sidechain from template: " << i << std::endl;
				} else {
					needToRepackAtAll = true;
				}
			} // for ( i in pose.size() )
		} // fa_input

		// create score function and add constraints for fullatom part
		if ( cmd_line_csts() ) {
			core::scoring::constraints::add_fa_constraints_from_cmdline(
				pose, *fa_scorefxn_
			);
		} else {
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(
				*fa_scorefxn_
			);
		}


		// Add coordinate constraints to non-loop regions if desired
		core::Real constrain_rigid_segments_weight = option[ OptionKeys::loops::constrain_rigid_segments ]();
		if ( constrain_rigid_segments_weight > 0.0 ) {
			protocols::loops::Loops coordconstraint_segments;
			//core::pose::Pose coordconstrainted_pose = pose;
			core::pose::Pose constraint_target_pose = pose;


			coordconstraint_segments = loops->invert( pose.size() );  // Invert the loops selection - i.e. the rigid areas are now defined
			std::cout << "Restraining the following segments: " << std::endl << coordconstraint_segments << std::endl;

			// ResidueTypeSet
			using namespace core;
			using namespace conformation;
			using namespace core::scoring;
			using namespace core::scoring::constraints;
			using namespace id;

			if ( pose.residue( pose.size() ).aa() != core::chemical::aa_vrt ) {
				pose.append_residue_by_jump(
					*ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ), pose.size()/2 );
			}

			//fpd  nmonomerres is #residues in a single subunit (excluding virtuals)
			core::Size rootres = pose.fold_tree().root();
			core::Size nmonomerres = pose.size()-1;
			core::conformation::symmetry::SymmetryInfoCOP symm_info;
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				core::conformation::symmetry::SymmetricConformation & SymmConf (
					dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
				symm_info = SymmConf.Symmetry_Info();
				nmonomerres = symm_info->num_independent_residues();
			}

			if ( !option[ OptionKeys::relax::coord_cst_width ].user() ) {
				Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
				// default is 0.5 (from idealize) -- maybe too small
				for ( Size i = 1; i<=nmonomerres; ++i ) {
					if ( !pose.residue(i).is_polymer() ) continue;
					if ( coordconstraint_segments.is_loop_residue( i ) ) {
						Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
						for ( Size ii = 1; ii<=nat_i_rsd.last_backbone_atom(); ++ii ) {
							pose.add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,rootres), nat_i_rsd.xyz( ii ),
								new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) ) );
						}

						// now cst symmetry mates
						// if (symm_info) {
						//  for ( core::conformation::symmetry::SymmetryInfo::Clones::const_iterator pos=symm_info->bb_clones( i ).begin(),
						//     epos=symm_info->bb_clones( i ).end(); pos != epos; ++pos ) {
						//   for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						//    pose.add_constraint( new CoordinateConstraint( AtomID(ii,*pos), AtomID(1,nres), nat_i_rsd.xyz( ii ),
						//          new HarmonicFunc( 0.0, coord_sdev ) ) );
						//   }
						//  }
						// }
					}
				}
			} else {
				Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
				Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ]() );
				for ( Size i = 1; i<=nmonomerres; ++i ) {
					if ( !pose.residue(i).is_polymer() ) continue;
					if ( coordconstraint_segments.is_loop_residue( i ) ) {
						Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
						for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
							pose.add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,rootres), nat_i_rsd.xyz( ii ),
								new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
						}
						// now cst symmetry mates
						// if (symm_info) {
						//  for ( core::conformation::symmetry::SymmetryInfo::Clones::const_iterator pos=symm_info->bb_clones( i ).begin(),
						//     epos=symm_info->bb_clones( i ).end(); pos != epos; ++pos ) {
						//   for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						//    pose.add_constraint( new CoordinateConstraint( AtomID(ii,*pos), AtomID(1,nres), nat_i_rsd.xyz( ii ),
						//          new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
						//   }
						//  }
						// }
					}
				}
			}


			fa_scorefxn_->set_weight( coordinate_constraint, constrain_rigid_segments_weight );
		}


		// ----------------------------------------------------------


		// do the same (again) for fit-to-density
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *fa_scorefxn_ );
		}

		if ( debug ) pose.dump_pdb(curr_job_tag + "_before_repack.pdb");
		if ( needToRepackAtAll ) { // kic refine does its own initial repacking
			TR << "Repacking required" << std::endl;
			TR << "Detecting disulfides" << std::endl;
			TR << "Annotated sequence before repack: " << pose.annotated_sequence(true) << std::endl;
			// repack loop + missing-density residues
			pose.conformation().detect_disulfides();
			using namespace core::pack::task;
			using namespace core::pack::task::operation;
			TaskFactoryOP tf = new TaskFactory;
			tf->push_back( new NoRepackDisulfides );
			tf->push_back( new InitializeFromCommandline );
			tf->push_back( new IncludeCurrent );
			tf->push_back( new RestrictToRepacking );
			PackerTaskOP taskstd = tf->create_task_and_apply_taskoperations( pose );
			core::pose::symmetry::make_residue_mask_symmetric( pose, needToRepack );
			// does nothing if pose is not symm
			taskstd->restrict_to_residues(needToRepack);

			fa_scorefxn_->show_line( TR, pose );
			core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );  //fpd symmetrize this
			mm->set_bb( false );
			mm->set_chi( true );
			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				if ( pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
					TR<<"disabling minimization on disulfide residue "<<i<<std::endl;
					mm->set_chi( i, false );
				}
			}

			//fpd symmetrize this
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				minimization_packing::symmetry::SymPackRotamersMover pack1( fa_scorefxn_, taskstd );
				pack1.apply( pose );

				core::optimization::symmetry::SymAtomTreeMinimizer mzr;
				core::optimization::MinimizerOptions options("lbfgs_armijo_nonmonotone", 1e-5, true, false);
				core::pose::symmetry::make_symmetric_movemap( pose, *mm );

				mzr.run( pose, *mm, *fa_scorefxn_, options );
			} else {
				protocols::minimization_packing::PackRotamersMover pack1( fa_scorefxn_, taskstd );
				pack1.apply( pose );

				// quick SC minimization
				core::optimization::AtomTreeMinimizer mzr;
				core::optimization::MinimizerOptions options("lbfgs_armijo_nonmonotone", 1e-5, true, false);
				mzr.run( pose, *mm, *fa_scorefxn_, options );
			}

			fa_scorefxn_->show_line( TR, pose );
		} else {
			//fpd
			(*fa_scorefxn_)(pose);
			TR << "No repacking required" << std::endl;
			(*fa_scorefxn_)(pose);
		}

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_repack.pdb");
	} // fullatom_output

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  intermediate relax the structure
	////
	////
	if ( intermedrelax() != "no" && (all_loops_closed) ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Intermediate Relax  " << std::endl;
		TR << "===" << std::endl;

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( option[ OptionKeys::loops::relax_with_foldtree ].user() ) {
			loops::fold_tree_from_loops( pose, *loops, f_new );
			pose.fold_tree( f_new );
			loops::add_cutpoint_variants( pose );
			fa_scorefxn_->set_weight( core::scoring::chainbreak,        option[ OptionKeys::relax::chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, option[ OptionKeys::relax::linear_chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, option[ OptionKeys::relax::overlap_chainbreak_weight ]() );
		}


		TR << pose.fold_tree() << std::endl;
		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
				setPoseExtraScore( pose, "brlx_loopcarms",  loops::loop_rmsd(native_pose_super, pose, *loops,true ) );
			}
		}

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "relax", true , true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_relax.pdb");
			if ( intermedrelax() == "relax" ) {
				relax::RelaxProtocolBaseOP relax_prot = relax::generate_relax_from_cmd();
				relax_prot->set_current_tag( curr_job_tag );
				relax_prot->apply( pose );
			} else if ( ( intermedrelax() == "fastrelax" ) ||
					( intermedrelax() == "seqrelax" ) ) {
				relax::FastRelax seqrelax( fa_scorefxn_, option[ OptionKeys::relax::sequence_file ]() );
				seqrelax.set_current_tag( curr_job_tag );
				seqrelax.apply( pose );
			}
			checkpoints_.checkpoint( pose, curr_job_tag, "relax", true);
		}
		checkpoints_.debug( curr_job_tag, "relax", (*fa_scorefxn_)( pose ) );
	} // intermediate relax the structure

	if ( intermedrelax() != "no" && (!all_loops_closed) ) {
		//The following keeps the score lines equivalent in the silent file.
		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
				setPoseExtraScore( pose, "brlx_loopcarms",  loops::loop_rmsd(native_pose_super, pose, *loops,true ) );
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////
	////
	////  Loop refine (fullatom type loop modelling )
	////
	////

	if ( refine() != "no" && (all_loops_closed) ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Refine " << std::endl;
		TR << "===" << std::endl;

		long starttime = time(NULL);

		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "bref_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "bref_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "bref_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "bref_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
			}
		}

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( refine() == "refine_kic" ) {
			loops::fold_tree_from_loops( pose, *loops, f_new, true /* include terminal cutpoints */);
		} else {
			loops::fold_tree_from_loops( pose, *loops, f_new);
		}
		pose.fold_tree( f_new );
		TR << "fold_tree_before_refine " << pose.fold_tree() << std::endl;

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "refine", true , true) ) {

			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_refine.pdb");
			if ( refine() == "refine_ccd" ) {
				loops::loop_mover::refine::LoopMover_Refine_CCD refine_ccd( loops, fa_scorefxn_ );
				refine_ccd.set_native_pose( new core::pose::Pose ( native_pose ) );
				refine_ccd.apply( pose );
			} else if ( refine() == "refine_kic" ) {
				//loops.remove_terminal_loops( pose );
				loops::loop_mover::refine::LoopMover_Refine_KIC refine_kic( loops, fa_scorefxn_ );
				refine_kic.set_native_pose( new core::pose::Pose ( native_pose ) );
				refine_kic.apply( pose );
			} else if ( refine() == "refine_kic_refactor" ) {
				using protocols::loop_modeling::LoopProtocol;
				using protocols::loop_modeling::LoopProtocolOP;
				using protocols::loop_modeling::LoopMover;
				using protocols::loop_modeling::LoopMoverOP;
				using protocols::loop_modeling::samplers::KicSampler;
				using protocols::loop_modeling::refiners::RepackingRefiner;
				using protocols::loop_modeling::refiners::RotamerTrialsRefiner;
				using protocols::loop_modeling::refiners::LocalMinimizationRefiner;
				using protocols::loop_modeling::utilities::PeriodicTask;
				using protocols::loop_modeling::loggers::LoggerOP;
				using protocols::loop_modeling::loggers::ProgressBar;
				using protocols::loop_modeling::loggers::ScoreVsRmsd;
				using protocols::loop_modeling::loggers::AcceptanceRates;

				LoopProtocolOP protocol = new LoopProtocol;
				LoopMoverOP mover = new LoopMover;
				LoggerOP acceptance_logger = new AcceptanceRates;
				loops::Loop loop = (*loops)[1];

				Size repack_period = 20;
				if ( option[OptionKeys::loops::repack_period].user() ) {
					repack_period = option[OptionKeys::loops::repack_period]();
				}

				mover->add_task(new KicSampler(acceptance_logger));
				mover->add_task(new PeriodicTask(new RepackingRefiner, repack_period));
				mover->add_task(new RotamerTrialsRefiner);
				mover->add_task(new LocalMinimizationRefiner);

				Size outer_cycles = 3;
				Size inner_cycles = 10 * loop.size();

				if ( option[OptionKeys::loops::outer_cycles].user() ) {
					outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
				}
				if ( option[OptionKeys::loops::max_inner_cycles].user() ) {
					Size max_cycles = option[OptionKeys::loops::max_inner_cycles]();
					inner_cycles = std::max(inner_cycles, max_cycles);
				}
				if ( option[OptionKeys::loops::fast] ) {
					outer_cycles = 3;
					inner_cycles = 12;
				}

				protocol->set_mover(mover);
				protocol->set_loop(loop);
				protocol->set_score_function(fa_scorefxn_);
				protocol->set_iterations(outer_cycles, inner_cycles, 2);
				protocol->add_logger(new ProgressBar);
				protocol->add_logger(new ScoreVsRmsd(native_pose, loop));
				protocol->add_logger(acceptance_logger);

				protocol->apply(pose);
			}

			if ( debug ) pose.dump_pdb(curr_job_tag + "_after_refine.pdb");
			checkpoints_.checkpoint( pose, curr_job_tag, "refine", true);

			// need to get the chainbreak score before the cutpoint variants are removed

			if ( option[ OptionKeys::loops::kic_use_linear_chainbreak ]() ) {
				setPoseExtraScore(
					pose, "final_chainbreak",
					pose.energies().total_energies()[ core::scoring::linear_chainbreak ]
				);
			} else {
				setPoseExtraScore(
					pose, "final_chainbreak",
					pose.energies().total_energies()[ core::scoring::chainbreak ]
				);
			}
		}

		checkpoints_.debug( curr_job_tag, "refine", (*fa_scorefxn_)( pose ) );

		loops::remove_cutpoint_variants( pose );
		// restore simple fold tree
		pose.fold_tree( f_orig );

		endtime = time(NULL);

		TR << "Refinetime: " << endtime - starttime << std::endl;

	}
	if ( refine() != "no" && (!all_loops_closed) ) {
		//The following keeps the score lines equivalent in the silent file.
		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
			}
		}
	}

	// if ( refine != "no" )

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  Maybe idealize the structure before relax ?
	////
	////
	if ( option[ OptionKeys::loops::idealize_after_loop_close ].user() && (all_loops_closed) ) {

		if ( debug ) {
			pose.dump_pdb(curr_job_tag + "_before_idealize.pdb");
		}
		protocols::idealize::IdealizeMover idealizer;
		idealizer.fast( false );
		idealizer.apply( pose );
		if ( debug ) {
			pose.dump_pdb(curr_job_tag + "_after_idealize.pdb");
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////
	////
	////  Relax the structure
	////
	////

	if ( relax() != "no" && (all_loops_closed) ) {
		TR << "====================================================================================" << std::endl;
		TR << "===" << std::endl;
		TR << "===   Relax  " << std::endl;
		TR << "===" << std::endl;

		core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
		if ( option[ OptionKeys::loops::relax_with_foldtree ].user() ) {
			loops::fold_tree_from_loops( pose, *loops, f_new );
			pose.fold_tree( f_new );
			loops::add_cutpoint_variants( pose );
			fa_scorefxn_->set_weight( core::scoring::chainbreak,        option[ OptionKeys::relax::chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, option[ OptionKeys::relax::linear_chainbreak_weight ]() );
			fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, option[ OptionKeys::relax::overlap_chainbreak_weight ]() );
		}


		TR << pose.fold_tree() << std::endl;
		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
			}
		}

		if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "relax", true , true) ) {
			if ( debug ) pose.dump_pdb(curr_job_tag + "_before_relax.pdb");
			if ( relax() == "relax" ) {
				relax::RelaxProtocolBaseOP relax_prot = relax::generate_relax_from_cmd();
				relax_prot->set_current_tag( curr_job_tag );
				relax_prot->apply( pose );
			} else if ( ( relax() == "fastrelax" ) ||
					( relax() == "seqrelax" ) ) {
				relax::FastRelax seqrelax( fa_scorefxn_, option[ OptionKeys::relax::sequence_file ]() );
				seqrelax.set_current_tag( curr_job_tag );
				seqrelax.apply( pose );
			} else if ( relax() == "minirelax" ) {
				protocols::relax::MiniRelax mini_relax( fa_scorefxn_ );
				mini_relax.set_current_tag( curr_job_tag );
				mini_relax.apply( pose );
			}
			checkpoints_.checkpoint( pose, curr_job_tag, "relax", true);
		}
		checkpoints_.debug( curr_job_tag, "relax", (*fa_scorefxn_)( pose ) );

		// restore simple fold tree
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_orig );

		fa_scorefxn_->set_weight( core::scoring::chainbreak, 0.0 );
		fa_scorefxn_->set_weight( core::scoring::linear_chainbreak, 0.0 );
		fa_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 0.0 );

		if ( debug ) pose.dump_pdb(curr_job_tag + "_after_relax.pdb");

		if ( option[ OptionKeys::loops::final_clean_fastrelax ]() ) {
			fa_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::constant_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::coordinate_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::angle_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
			fa_scorefxn_->set_weight( core::scoring::big_bin_constraint, 0.0 );
			if ( !checkpoints_.recover_checkpoint( pose, curr_job_tag, "ffrelax", true , true ) ) {
				relax::FastRelax fast_relax( fa_scorefxn_ );
				fast_relax.set_current_tag( curr_job_tag );
				fast_relax.apply( pose );
				checkpoints_.checkpoint( pose, curr_job_tag, "ffrelax", true);
			}
			checkpoints_.debug( curr_job_tag, "ffrelax", (*fa_scorefxn_)( pose ) );
		}
	} // relax the structure
	if ( relax() != "no" && (!all_loops_closed) ) {
		//The following keeps the score lines equivalent in the silent file.
		if ( compute_rmsd() ) {
			setPoseExtraScore( pose, "brlx_irms",  core::scoring::CA_rmsd( start_pose, pose ) );
			if ( option[ in::file::native ].user() ) {
				if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
					core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
				}
				setPoseExtraScore( pose, "brlx_rms",   core::scoring::native_CA_rmsd(native_pose, pose ) );
				setPoseExtraScore( pose, "brlx_corerms", native_loop_core_CA_rmsd(native_pose, pose, *loops, corelength ) );
				setPoseExtraScore( pose, "brlx_looprms",  loops::loop_rmsd(native_pose_super, pose, *loops ) );
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
	////
	////   Statistics
	////

	TR << "====================================================================================" << std::endl;
	TR << "===" << std::endl;
	TR << "===  Getting Statistics " << std::endl;
	TR << "===" << std::endl;
	TR << "===" << std::endl;

	if ( compute_rmsd() ) {
		setPoseExtraScore(
			pose, "irms",  core::scoring::CA_rmsd( start_pose, pose )
		);
		if ( option[ in::file::native ].user() ) {
			if (  option[ OptionKeys::loops::superimpose_native ]()  ) {
				core::scoring::superimpose_pose( native_pose_super, pose, atom_map );
			}
			setPoseExtraScore( pose, "rms",     core::scoring::native_CA_rmsd( native_pose, pose ));
			setPoseExtraScore( pose, "looprms", loops::loop_rmsd( native_pose_super, pose, *loops, false /*CA_only*/, true /*bb_only*/ ));
			setPoseExtraScore( pose, "loop_heavy_rms", loops::loop_rmsd( native_pose_super, pose, *loops, false /*CA_only*/, false /*bb_only*/ ));
			setPoseExtraScore( pose, "loopcarms", loops::loop_rmsd( native_pose_super, pose, *loops, true ));
			setPoseExtraScore( pose, "corerms", native_loop_core_CA_rmsd( native_pose, pose, *loops, corelength ) );
			setPoseExtraScore( pose, "corelen", corelength );
			//if( pose.is_fullatom() ) addScoresForLoopParts( pose, loops, (*fa_scorefxn_), native_pose, 10 );
			//else                     addScoresForLoopParts( pose, loops, (*cen_scorefxn_), native_pose, 10 );
		}
	}

	core::Real final_score;

	if ( debug ) pose.dump_pdb(curr_job_tag + "_before_final_rescore.pdb");

	if ( fullatom_output ) final_score = (*fa_scorefxn_)(pose);  // may include constraint score
	else                 final_score = (*cen_scorefxn_)(pose); // may include constraint score


	core::pose::setPoseExtraScore(
		pose, std::string("final_looprelax_score"), final_score
	);

} // LoopRelaxMover

// {{{1
/// @brief Must be called before the Loops data can be read from.
void LoopRelaxMover::resolve_loopfile_indices( core::pose::Pose const & pose )
{
	guarded_loops_->resolve_loop_indices( pose );
}

// {{{1
/// @details By setting the loops object directly, the requirement is relaxed that a Pose be first
/// given to the LoopRelax object (through a call to apply() or resolve_loopfile_indices)
/// before the get_loops() function may be called.
void LoopRelaxMover::loops( protocols::loops::LoopsOP const val ) {
	//there is an assertion assert( !in_charge_ ) in set_loops_pointer()
	// as the guraded_loops_ is made with in_charge_ = true in the Constructor
	// and this function is called in parse_my_tags to set the loops, the assertion
	// seems to be wrong, or the use of this class is wrong.
	// I have no idea what this assertion is supposed to achieve...
	guarded_loops_->in_charge( false );
	guarded_loops_->set_loops_pointer( val );
	guarded_loops_->in_charge( true );
}

// {{{1
/// @brief Set the loop file data.  This will require that
void LoopRelaxMover::loops_file_data( loops::LoopsFileData const & loopfiledata ) {
	guarded_loops_->loops( loopfiledata );
}

// {{{1
protocols::loops::LoopsCOP
LoopRelaxMover::get_loops() const {
	return guarded_loops_->loops();
}

// {{{1
protocols::loops::LoopsOP
LoopRelaxMover::get_loops() {
	return guarded_loops_->loops();
}


// {{{1
std::string
LoopRelaxMover::get_name() const {
	return "LoopRelaxMover";
}

// {{{1
void LoopRelaxMover::set_defaults_() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	cmd_line_csts_   = true;
	copy_sidechains_ = true;
	remodel_         = option[ OptionKeys::loops::remodel ]();
	intermedrelax_   = option[ OptionKeys::loops::intermedrelax ]();
	refine_          = option[ OptionKeys::loops::refine ]() ;
	relax_           = option[ OptionKeys::loops::relax ]();
	n_rebuild_tries_ = 3;
	rebuild_filter_  = 999.0;
	compute_rmsd( true );

	// use score4L by default (will be symm if needed)
	cen_scorefxn_ = loops::get_cen_scorefxn();
	loops::loop_mover::loops_set_chainbreak_weight(  cen_scorefxn_, 1 );

	// get cmd line scorefxn by default (will be symm if needed)
	fa_scorefxn_ = loops::get_fa_scorefxn();
}

// {{{1
void LoopRelaxMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data )
{
	using namespace utility::tag;
	using namespace protocols::moves;
	cmd_line_csts( tag->getOption< bool >( "cmd_line_csts", true ) );
	copy_sidechains( tag->getOption< bool >( "copy_sidechains", true ) );
	rebuild_filter( tag->getOption< core::Real >( "rebuild_filter", 999 ) );
	compute_rmsd( tag->getOption< bool >( "compute_rmsd", true ) );
	n_rebuild_tries( tag->getOption< core::Size >( "n_rebuild_tries", 10 ) );
	remodel( tag->getOption< std::string >( "remodel", "quick_ccd" ) );
	intermedrelax( tag->getOption< std::string >( "intermedrelax", "no" ) );
	refine( tag->getOption< std::string >( "refine", "refine_ccd" ) );
	relax( tag->getOption< std::string >( "relax", "no" ) );
	if ( tag->getOption< bool >( "read_fragments", false ) ) {
		loops::read_loop_fragments( frag_libs_ );
	}

	//read tag "scorefxn"
	fa_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	cen_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data, "score4L" ) );
	//cen_scorefxn_->set_weight( core::scoring::chainbreak, 10.0 / 3.0 );
	loops::loop_mover::loops_set_chainbreak_weight(  cen_scorefxn_, 1 );

	loops( loops::loops_definers::load_loop_definitions(tag, data, pose) );
}

// {{{1
std::string LoopRelaxMoverCreator::keyname() const
{
	return LoopRelaxMoverCreator::mover_name();
}

// {{{1
protocols::moves::MoverOP LoopRelaxMoverCreator::create_mover() const {
	return new LoopRelaxMover;
}

// {{{1
std::string LoopRelaxMoverCreator::mover_name()
{
	return "LoopRelaxMover";
}

// {{{1
protocols::moves::MoverOP LoopRelaxMover::fresh_instance() const {
	return utility::pointer::make_shared< LoopRelaxMover >();
}

// {{{1
protocols::moves::MoverOP LoopRelaxMover::clone() const {
	return utility::pointer::make_shared< LoopRelaxMover >( *this );
}
// }}}1

} // namespace loops
} // namespace protocols
