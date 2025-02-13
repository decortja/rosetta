// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/RealMetric.cc
///
/// @brief Main class for simple metrics.
/// @author Jared Adolf-Bryfogle ( jadolfbr@gmail.com )

// Unit Headers
#include <core/simple_metrics/RealMetric.hh>

// Protocol Headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/util.hh>
#include <core/simple_metrics/SimpleMetricData.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {

RealMetric::RealMetric():
	SimpleMetric("RealMetric")
{}


RealMetric::~RealMetric() = default;

RealMetric::RealMetric( RealMetric const & src ):
	SimpleMetric( src )
{}

RealMetric &
RealMetric::operator=( RealMetric const & ot ) {
	SimpleMetric::operator=( ot );
	return *this;
}

void
RealMetric::apply( std::string const & out_tag, pose::Pose & pose, bool override_existing ) const {

	core::Real value = calculate( pose );
	MetricKey mk;
	core::Real stored_value;

	if ( ( ! override_existing ) && get_sm_data(pose)->get_value(out_tag, stored_value) ) {
		throw_sm_override_error(out_tag, name());
	}
	get_sm_data(pose)->set_value(mk, out_tag, value);
}

core::Real
RealMetric::cached_calculate(pose::Pose const & pose, bool use_cache, std::string prefix, std::string suffix, bool fail_on_missing_cache) const {
	std::string name = prefix + get_final_sm_type() + suffix;

	if ( use_cache && has_sm_data( pose ) ) {
		core::Real value;
		bool data_found = get_sm_data(pose)->get_value(name, value);
		if ( data_found ) {
			return value;
		} else if ( fail_on_missing_cache ) {
			utility_exit_with_message("Could not find RealMetric: "+name+" in pose");
		} else {
			return calculate(pose);
		}
	} else {
		return calculate(pose);
	}
}

utility::vector1< std::string >
RealMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

} //namespace simple_metrics
} //namespace core

#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::RealMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric>( this ) );

}

template< class Archive >
void
core::simple_metrics::RealMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric >( this ) );


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::RealMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::RealMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_RealMetric )
#endif // SERIALIZATION

