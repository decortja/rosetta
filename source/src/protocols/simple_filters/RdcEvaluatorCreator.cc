// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/RdcEvaluatorCreator.hh
/// @brief  Header for RdcEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RdcEvaluatorCreator.hh>

// Package Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RDC_Evaluator.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//// C++ headers

// due to template function


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>


//Auto Headers


#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static basic::Tracer tr( "protocols.evalution.RdcEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RdcEvaluatorCreator::~RdcEvaluatorCreator() = default;

void RdcEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::rdc );

}

void RdcEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::rdc ].user() ) {
		using RdcVector = utility::vector1<std::string>;
		RdcVector const& rdc( option[ OptionKeys::evaluation::rdc ]() );
		utility::vector1< core::Size> empty_selection;
		for ( auto it=rdc.begin(); it!=rdc.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != rdc.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
					"need to specify dupletts <rdcs> <column> with option -evaluation:rdc   last read: "+fname );
			}
			eval.add_evaluation( utility::pointer::make_shared< simple_filters::SelectRDC_Evaluator >( empty_selection, column, fname ) );
		} // iterate over tripletts in option -rmsd
	}
	if ( option[ OptionKeys::evaluation::built_in_rdc ].user() ) {
		eval.add_evaluation( utility::pointer::make_shared< simple_filters::RDC_Evaluator >(option[ OptionKeys::evaluation::built_in_rdc ]()) );
	} // iterate over tripletts in option -rmsd
}

std::string RdcEvaluatorCreator::type_name() const {
	return "RdcEvaluatorCreator";
}

} //namespace
} //namespace
