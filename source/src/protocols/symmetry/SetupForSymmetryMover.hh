// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_symmetry_SetupForSymmetryMover_hh
#define INCLUDED_protocols_simple_moves_symmetry_SetupForSymmetryMover_hh

// Unit headers
#include <protocols/symmetry/SetupForSymmetryMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>

#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

// Utility Headers

namespace protocols {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////

class SetupForSymmetryMover : public protocols::moves::Mover
{
public:

	/// @brief default constructor, reads from the global options object
	SetupForSymmetryMover();

	/// @brief reads from a (possibly-local) option collection object
	SetupForSymmetryMover( utility::options::OptionCollection const & options );

	SetupForSymmetryMover( core::conformation::symmetry::SymmDataOP symmdata );

	SetupForSymmetryMover(
		core::conformation::symmetry::SymmDataOP symmdata,
		utility::options::OptionCollection const & options
	);

	SetupForSymmetryMover( std::string const & symmdef_file );
	SetupForSymmetryMover( std::string const & symmdef_file,
		utility::options::OptionCollection const & options );

	~SetupForSymmetryMover() override;

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< SetupForSymmetryMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const override;

	// setter
	void slide_into_contact(bool val) { slide_ = val; }

	/// @brief Sets whether or not the input asymmetric pose's datacache should be copied into
	///        the new symmetric pose.
	/// @param[in] preserve_cache If true, input pose's datacache is copied into new symmetric pose
	///                           If false, input pose's datacache is cleared (default = false)
	void
	set_preserve_datacache( bool const preserve_cache );

	void
	set_keep_pdb_info_labels( bool const keep_pdb_info_labels );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void process_symmdef_file(std::string tag);

	static
	void
	options_read_in_ctor( utility::options::OptionKeyList & opts );

	void
	set_refinable_lattice( bool setting );

private:

	/// @brief   constructs a symmetric pose with a symmetric conformation and energies object.
	/// @details Calls core::pose::make_symmetric_pose().  If preserve_datacache is set, this
	///          also copies the datacache into the new symmetric pose.
	/// @param[in,out] pose Input asymmetric pose, output symmetric pose
	void
	make_symmetric_pose( core::pose::Pose & pose ) const;

	void
	read_refinable_lattice( utility::options::OptionCollection const & options );

private:
	bool slide_;
	bool cryst1_; //use cryst1 line
	bool preserve_datacache_;
	core::conformation::symmetry::SymmDataOP symmdef_;
	std::string symdef_fname_from_options_system_;
	bool refinable_lattice_was_set_;
	bool refinable_lattice_;
	bool keep_pdb_info_labels_;
	bool set_global_symmetry_at_parsetime_; //Flag to set if the symdef option should be set globally.
};

///////////////

///@brief Extract the Asymmetric unit only from the pose.  This is the pose that would have existed before symmetrization.
class ExtractAsymmetricUnitMover : public protocols::moves::Mover
{
public:

	// default constructor
	ExtractAsymmetricUnitMover();

	~ExtractAsymmetricUnitMover() override;

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< ExtractAsymmetricUnitMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;


public:
	/// @brief if keep_virtual_residues is true, virtual residues will remain in the pose, otherwise
	///        they will be removed
	void
	set_keep_virtual_residues( bool const keep_virt );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const override;


private:
	/// @brief if keep_virtual_residues is true, virtual residues will remain in the pose, otherwise
	///        they will be removed
	bool keep_virtual_residues_;

};

class ExtractAsymmetricPoseMover : public protocols::moves::Mover
{
public:

	// default constructor
	ExtractAsymmetricPoseMover();

	~ExtractAsymmetricPoseMover() override;

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< ExtractAsymmetricPoseMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;


public:
	/// @brief if true, clears symmetry_definition option from pose.
	void
	clear_sym_def( bool const clear_sym_def );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const override;

private:
	/// @brief if true, clears symmetry_definition option from pose.
	bool clear_sym_def_;
};


}
} // rosetta
#endif
