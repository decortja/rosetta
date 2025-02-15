// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/AtomPairConstraintGenerator.cc
/// @brief Generates atom pair constraints for a set of residues from the current or reference pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/AtomPairConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
static basic::Tracer TR( "protocols.constraint_generator.AtomPairConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
AtomPairConstraintGeneratorCreator::create_constraint_generator() const
{
	return utility::pointer::make_shared< AtomPairConstraintGenerator >();
}

std::string
AtomPairConstraintGeneratorCreator::keyname() const
{
	return AtomPairConstraintGenerator::class_name();
}

AtomPairConstraintGenerator::AtomPairConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( AtomPairConstraintGenerator::class_name() ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	reference_pose_(),
	sd_( 0.5 ),
	weight_( 1.0 ),
	ca_only_( true ),
	use_harmonic_( false ),
	unweighted_( false ),
	max_distance_( 12.0 ),
	min_seq_sep_( 8 )
{
	core::select::residue_selector::ResidueSelectorCOP true_selector( new core::select::residue_selector::TrueResidueSelector );
	core::select::residue_selector::NotResidueSelector false_selector( true_selector );
	set_secondary_residue_selector( false_selector );
}

AtomPairConstraintGenerator::~AtomPairConstraintGenerator() = default;

protocols::constraint_generator::ConstraintGeneratorOP
AtomPairConstraintGenerator::clone() const
{
	return utility::pointer::make_shared< AtomPairConstraintGenerator >( *this );
}

void
AtomPairConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	bool const use_native = tag->getOption< bool >( "native", false );
	if ( use_native ) {
		set_reference_pose( get_native_pose() );
		if ( ! reference_pose_ ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "'native' option for AtomPairConstraintGenerator specified, but no native pose is availible." );
		}
	}

	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( *selector );

	std::string const selector_name = tag->getOption< std::string >( "secondary_selector", "" );
	if ( !selector_name.empty() ) secondary_selector_ = core::select::residue_selector::get_residue_selector( selector_name, data );

	set_sd( tag->getOption< core::Real >( "sd", sd_ ) );
	set_weight( tag->getOption< core::Real >( "weight", weight_ ) );
	set_ca_only( tag->getOption< bool >( "ca_only", ca_only_ ) );
	set_max_distance( tag->getOption< core::Real >( "max_distance", max_distance_ ) );
	set_min_seq_sep( tag->getOption< core::Size >( "min_seq_sep", min_seq_sep_ ) );
	set_use_harmonic_function( tag->getOption< bool >( "use_harmonic", use_harmonic_ ) );
	set_unweighted_function( tag->getOption< bool >( "unweighted", unweighted_ ) );

	if ( !selector_ ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "AtomPairConstraintGenerator requires a residue selector, but one is not set.\n" );
	}
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	debug_assert( selector_ );
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	core::select::residue_selector::ResidueSubset const secondary_subset = secondary_selector_->apply( pose );
	return generate_atom_pair_constraints( pose, subset, secondary_subset );
}

/// @brief Provide citations to the passed CitationCollectionList.
/// This allows the constraint generator to provide citations for itself
/// and for any modules that it invokes.
/// @details Cites Tom Linsky, who created the constraint generator framework.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
AtomPairConstraintGenerator::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP citation(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		class_name(), CitedModuleType::ConstraintGenerator,
		"Thomas W. Linsky", "Neoleukin Therapeutics", "tlinsky@gmail.com",
		"Created the ConstraintGenerator framework and the AtomPairConstraintGenerator."
		)
	);
	citations.add(citation);
}

void
AtomPairConstraintGenerator::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}

void
AtomPairConstraintGenerator::set_secondary_residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	secondary_selector_ = selector.clone();
}

void
AtomPairConstraintGenerator::set_sd( core::Real const sd )
{
	sd_ = sd;
}

void
AtomPairConstraintGenerator::set_ca_only( bool const ca_only )
{
	ca_only_ = ca_only;
}

void
AtomPairConstraintGenerator::set_use_harmonic_function( bool const use_harmonic )
{
	use_harmonic_ = use_harmonic;
}

void
AtomPairConstraintGenerator::set_unweighted_function( bool const unweighted )
{
	unweighted_ = unweighted;
}

void
AtomPairConstraintGenerator::set_max_distance( core::Real const max_dist )
{
	max_distance_ = max_dist;
}

void
AtomPairConstraintGenerator::set_min_seq_sep( core::Size const min_seq_sep )
{
	min_seq_sep_ = min_seq_sep;
}

void
AtomPairConstraintGenerator::set_weight( core::Real const weight )
{
	weight_ = weight;
}

void
AtomPairConstraintGenerator::set_reference_pose( core::pose::PoseCOP ref_pose )
{
	reference_pose_ = ref_pose;
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::select::residue_selector::ResidueSubset const & subset2 ) const
{
	if ( reference_pose_ ) {
		return generate_atom_pair_constraints( pose, *reference_pose_, subset, subset2 );
	} else {
		return generate_atom_pair_constraints( pose, pose, subset, subset2 );
	}
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::select::residue_selector::ResidueSubset const & subset2 ) const
{
	core::id::SequenceMapping const seqmap = create_sequence_mapping( pose, ref_pose );
	if ( core::select::residue_selector::has_any_true_selection( subset2 ) ) {
		return generate_atom_pair_constraints( pose, ref_pose, subset, subset2, seqmap );
	}
	return generate_atom_pair_constraints( pose, ref_pose, subset, seqmap );
}

core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::select::residue_selector::ResidueSubset const & subset2,
	core::id::SequenceMapping const & seqmap ) const
{
	core::scoring::constraints::ConstraintCOPs csts;

	core::Size const nres = compute_nres_in_asymmetric_unit( pose );
	for ( core::Size ires=1; ires<=nres; ++ires ) {
		if ( !subset[ ires ] ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
		core::Size const ref_ires = seqmap[ ires ];
		if ( ref_ires == 0 ) {
			TR.Debug << "Residue " << ires << " not found in reference pose, skipping" << std::endl;
			continue;
		}
		MappedAtoms const iatoms = atoms_to_constrain( pose.residue( ires ), ref_pose.residue( ref_ires ) );

		for ( core::Size jres=1; jres<=pose.size(); ++jres ) {
			if ( !subset2[ jres ] ) continue;
			if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
			core::Size const ref_jres = seqmap[ jres ];
			if ( ref_jres == 0 ) {
				TR.Debug << "Residue " << jres << " not found in reference pose, skipping" << std::endl;
				continue;
			}

			MappedAtoms const jatoms = atoms_to_constrain( pose.residue( jres ), ref_pose.residue( ref_jres ) );
			add_constraints( csts, ires, jres, ref_pose.residue( ref_ires ), ref_pose.residue( ref_jres ), iatoms, jatoms );
		} // jres loop
	} // ires loop
	return csts;

}
core::scoring::constraints::ConstraintCOPs
AtomPairConstraintGenerator::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose,
	core::select::residue_selector::ResidueSubset const & subset,
	core::id::SequenceMapping const & seqmap ) const
{
	core::scoring::constraints::ConstraintCOPs csts;

	core::Size const nres = compute_nres_in_asymmetric_unit( pose );
	for ( core::Size ires=1; ires<=nres; ++ires ) {
		if ( !subset[ ires ] ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
		core::Size const ref_ires = seqmap[ ires ];
		if ( ref_ires == 0 ) {
			TR.Debug << "Residue " << ires << " not found in reference pose, skipping" << std::endl;
			continue;
		}
		MappedAtoms const iatoms = atoms_to_constrain( pose.residue( ires ), ref_pose.residue( ref_ires ) );

		for ( core::Size jres=ires+min_seq_sep_; jres<=pose.size(); ++jres ) {
			if ( !subset[ jres ] ) continue;
			if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
			core::Size const ref_jres = seqmap[ jres ];
			if ( ref_jres == 0 ) {
				TR.Debug << "Residue " << jres << " not found in reference pose, skipping" << std::endl;
				continue;
			}

			MappedAtoms const jatoms = atoms_to_constrain( pose.residue( jres ), ref_pose.residue( ref_jres ) );
			add_constraints( csts, ires, jres, ref_pose.residue( ref_ires ), ref_pose.residue( ref_jres ), iatoms, jatoms );
		} // jres loop
	} // ires loop
	return csts;
}

void
AtomPairConstraintGenerator::add_constraints(
	core::scoring::constraints::ConstraintCOPs & csts,
	core::Size const pose_resid1,
	core::Size const pose_resid2,
	core::conformation::Residue const & ref_ires,
	core::conformation::Residue const & ref_jres,
	MappedAtoms const & iatoms,
	MappedAtoms const & jatoms ) const
{
	for ( auto iatom : iatoms ) {
		for ( auto jatom : jatoms ) {
			core::Real const dist = ref_ires.xyz( iatom.ref_atom ).distance( ref_jres.xyz( jatom.ref_atom ) );
			if ( dist > max_distance_ ) continue;

			core::id::AtomID const pose_atom1( iatom.pose_atom, pose_resid1 );
			core::id::AtomID const pose_atom2( jatom.pose_atom, pose_resid2 );

			core::scoring::func::FuncOP pair_func;
			if ( !use_harmonic_ ) {  // This is the default behaviour
				pair_func = utility::pointer::make_shared< core::scoring::func::SOGFunc >( dist, sd_ );
			} else {
				pair_func = core::scoring::func::FuncOP (new core::scoring::func::HarmonicFunc( dist, sd_ ) );
			}

			core::scoring::constraints::ConstraintCOP newcst;
			if ( !unweighted_ ) {  // This is the default behaviour
				core::scoring::func::FuncOP weighted_func = scalar_weighted( pair_func, weight_ );
				newcst = utility::pointer::make_shared< core::scoring::constraints::AtomPairConstraint >( pose_atom1, pose_atom2, weighted_func );
			} else {
				newcst = utility::pointer::make_shared< core::scoring::constraints::AtomPairConstraint >( pose_atom1, pose_atom2, pair_func );
			}

			csts.push_back( newcst );
			TR.Debug << "atom_pair_constraint generated for residue " << pose_atom1 << " and " << pose_atom2 << ", distance=" << dist << " with weight " << weight_ << std::endl;
		}
	}
}

core::id::SequenceMapping
AtomPairConstraintGenerator::create_sequence_mapping( core::pose::Pose const & pose, core::pose::Pose const & ref_pose ) const
{
	bool const same_length = ( pose.size() == ref_pose.size() );
	bool const same_sequence = ( pose.sequence() == ref_pose.sequence() );

	if ( same_length && same_sequence ) {
		return core::id::SequenceMapping::identity( pose.size() );
	} else { // !same_sequence || !same_length
		TR << "Input structure and native differ in ";
		if ( !same_length ) TR << "length and sequence ";
		else if ( !same_sequence ) TR << "sequence ";
		TR << "- aligning on PDB identity or sequence." << std::endl;
		return core::pose::sequence_map_from_pdbinfo( pose, ref_pose );
	}
}

void
AtomPairConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "native", xsct_rosetta_bool, "Restrain to native distance?", "false" )
		+ XMLSchemaAttribute( "sd", xsct_real, "Standard deviation for distance constraint" )
		+ XMLSchemaAttribute( "weight", xsct_real, "Weight of distance constraint" )
		+ XMLSchemaAttribute( "ca_only", xsct_rosetta_bool, "Only make constraints between alpha carbons" )
		+ XMLSchemaAttribute::attribute_w_default( "use_harmonic", xsct_rosetta_bool, "If true, use harmonic function instead of SOG function", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "unweighted", xsct_rosetta_bool, "If true, SCALARWEIGHTEDFUNC is not added to the constraint definition", "false" )
		+ XMLSchemaAttribute( "max_distance", xsct_real, "Do not add constraints if atoms are farther apart than this" )
		+ XMLSchemaAttribute( "min_seq_sep", xsct_non_negative_integer, "Minimum sequence separation between constrained residues" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Selector specifying residues to be constrained. When not provided, all residues are selected" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "secondary_selector", "With a secondary selector, constraints are generated between the residues of primary selector vs. "
		"secondary selector. min_seq_seq does not apply here, but max_distance does." );
	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Generates atom pair constraints between specified residues",
		attlist );
}

void
AtomPairConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AtomPairConstraintGenerator::provide_xml_schema( xsd );
}




AtomPairConstraintGenerator::MappedAtoms
AtomPairConstraintGenerator::atoms_to_constrain(
	core::conformation::Residue const & pose_rsd,
	core::conformation::Residue const & ref_rsd ) const
{
	utility::vector1< MappedAtom > atoms;
	if ( pose_rsd.has( "CA" ) && ref_rsd.has( "CA" ) ) {
		atoms.push_back( MappedAtom( pose_rsd.type().atom_index( "CA" ), ref_rsd.type().atom_index( "CA" ) ) );
	}
	if ( !ca_only_ ) {
		if ( ( pose_rsd.nbr_atom() != atoms[1].pose_atom ) && ( ref_rsd.nbr_atom() != atoms[1].ref_atom ) ) {
			atoms.push_back( MappedAtom( pose_rsd.nbr_atom(), ref_rsd.nbr_atom() ) );
		}
	}
	return atoms;
}

} //protocols
} //constraint_generator
