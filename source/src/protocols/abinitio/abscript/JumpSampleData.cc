// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/JumpSampleData.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/JumpSampleData.hh>

// Package headers

// Project headers

// Utility Headers

// Basic Headers
#include <basic/Tracer.hh>

#include <protocols/abinitio/abscript/JumpSampleDataCreator.hh> // AUTO IWYU For JumpSampleDataCreator
#include <utility/excn/Exceptions.hh> // AUTO IWYU For BadInput, CREATE_EXCEPTION


//option includes

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr( "protocols.environment.movers.JumpSampleData", basic::t_info );

static std::string TYPE_NAME() { return "JumpSampleData"; }

namespace protocols {
namespace abinitio {
namespace abscript {

basic::datacache::WriteableCacheableDataOP
JumpSampleDataCreator::create_data( std::istream &in ) const {
	return utility::pointer::make_shared< JumpSampleData >( in );
}

std::string JumpSampleDataCreator::keyname() const{
	return TYPE_NAME();
}

JumpSampleData::JumpSampleData( std::istream &in ) :
	Parent()
{
	std::string token;

	in >> token;
	if ( token != "MOVERKEY" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "JumpSampleData tried to read an improperly formatted SilentFile remark." );
	}

	in >> moverkey_;

	core::kinematics::FoldTree ft;
	in >> ft;

	jump_sample_ = jumping::JumpSample( ft );
}

JumpSampleData::JumpSampleData( std::string const & moverkey,
	jumping::JumpSample const & jump_sample ) :
	Parent(),
	jump_sample_( jump_sample ),
	moverkey_( moverkey )
{}

void JumpSampleData::write( std::ostream &out ) const {

	out << TYPE_NAME() ;
	out << " MOVERKEY ";
	out << moverkey_ << " ";

	// get rid of nonsense new line in FT out method.
	std::stringstream ss;
	ss << jump_sample_.fold_tree();
	std::string ft( ss.str() );
	ft.erase( ft.length() - 1 );

	out << ft;

}

basic::datacache::CacheableDataOP
JumpSampleData::clone() const {
	return utility::pointer::make_shared< JumpSampleData >( *this );
}

std::string JumpSampleData::datatype() const {
	return TYPE_NAME();
}


} // abscript
} // abinitio
} // protocols
