// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/BluePrint.hh
/// @brief the basic idea of BluePrint is on the remodel Possu wrote in rosetta++.
/// @author Nobuyasu Koga (nobuyasu@uw.edu)

#ifndef INCLUDED_protocols_parser_BluePrint_HH
#define INCLUDED_protocols_parser_BluePrint_HH

// Unit  header
#include <protocols/parser/BluePrint.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>

#include <string>
#include <map>

#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace parser {


class BluePrint : public utility::VirtualBase {
public:


	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::MoveMapOP MoveMapOP;

public: // constructor/destructor


	/// @brief default constructor
	BluePrint();

	/// @brief value constructor
	BluePrint( String const & filename );

	/// @brief destructor
	~BluePrint() override;

	/// @brief copy constructor
	BluePrint( BluePrint const & rval );


public: // accessor


	/// @brief total residue number defined in blueprint file
	core::Size total_residue() const;

	/// @brief total residue number without ligand defined in blueprint file
	core::Size total_residue_wolig() const;

	/// @brief sequence defined in blueprint file
	String sequence() const;

	/// @brief a~mino acid type at a position in blueprint file
	char sequence( core::Size seqpos ) const;

	/// @brief secondary structures defined in blueprint file
	String secstruct() const;

	/// @brief secondary structure at a position in blueprint file
	char secstruct( core::Size seqpos ) const;

	/// @brief abego defined in bludprint file
	utility::vector1< String > abego() const;

	/// @brief secondary structure at a position in blueprint file
	String abego( core::Size seqpos ) const;

	/// @brief residue number at each position in blueprint file
	core::Size resnum( core::Size seqpos ) const;

	/// @brief translate residue number of pose to that of blueprint file
	core::Size resnum_map( core::Size resnum_pose ) const;

	/// @brief return build type at each position
	char buildtype( core::Size seqpos ) const;

	/// @brief return build type at each position
	String extra( core::Size seqpos ) const;

	String insertion( core::Size i ) const;

	/// @brief helix pairings defined at the line of HHPAIR in blueprint
	String helix_pairings() const;

	/// @brief strand pairings defined at the line of SSPAIR in blueprint
	String strand_pairings() const;

	/// @brief strand pairings defined at the line of SSPAIR in blueprint
	String hss_triplets() const;

	/// @brief secondary structure information
	//SS_Info2_OP ssinfo() const;


public: //


	/// @brief reading blueprint files
	bool read_blueprint( String const & );

	/// @brief read blueprint file from stream
	bool read_blueprint_stream( std::istream & data, std::string const & filename );

	/// @brief set secondary structure into pose
	void insert_ss_into_pose( Pose & pose );

	/// @brief set movemap based on blueprint
	void set_movemap( MoveMapOP & movemap );


private: // to be removed


	/// @brief set strand pairings
	/// removed StrandPairings set_strand_pairings( SS_Info2_OP const & ssinfo, StrandPairings const & spairs ) const;


private: // data

	/// @brief total residue number defined in blueprint
	core::Size total_residue_;

	/// @brief total residue without ligand
	core::Size total_residue_wolig_;

	/// @brief sequence defined in blueprint
	String sequence_;

	/// @brief amino acid sequence defined in blueprint
	String secstruct_;

	/// @brief residue number of each position in blueprint
	utility::vector1< core::Size > resnum_;

	/// @brief amino acid type at a position in blueprint
	utility::vector1< char > resname_;

	/// @brief secondary structure type defined in blueprint
	utility::vector1< char > sstype_;

	/// @brief abego type defined in blueprint
	utility::vector1< String > abego_;

	/// @brief build type at each position in blueprint
	utility::vector1< char > buildtype_;

	/// @brief extra infomation at each position in blueprint
	utility::vector1< String > extra_;

	/// @brief pdb file name for insertion
	utility::vector1< String > insertion_;

	/// @brief translate pose residue number to blueprint residue number
	std::map< core::Size, core::Size > resnum_map_;

	/// @brief secondary structure information
	// SS_Info2_OP ss_info_;

	/// @brief strand pairings defined at the line of SSPAIR in blueprint
	String strand_pairings_;

	/// @brief helix pairings defined at the line of HHPAIR in blueprint
	String helix_pairings_;

	/// @brief helix-strand-strand triple defined at the line of HSSTRIPLE in blueprint
	String hss_triplets_;


}; //BluePrint

} // parser
} // protocols

#endif
