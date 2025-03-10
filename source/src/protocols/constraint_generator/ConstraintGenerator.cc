// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGenerator.cc
///
/// @brief Abstract class for generating constraints
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit Headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Protocol Headers

// Core headers

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.moves.ConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

ConstraintGenerator::ConstraintGenerator( std::string const & class_name ):
	utility::VirtualBase(),
	class_name_( class_name ),
	id_( "unnamed_constraint_generator" )
{}

ConstraintGenerator::~ConstraintGenerator() = default;

/// @brief Provide citations to the passed CitationCollectionList.
/// This allows the constraint generator to provide citations for itself
/// and for any modules that it invokes.
/// @details This base class version does nothing.  It should be overridden by derived classes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
ConstraintGenerator::provide_citation_info(
	basic::citation_manager::CitationCollectionList &
) const {
	//GNDN.
}

/// @brief This is called if this mover is instantiated from XML
void
ConstraintGenerator::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	// if there are any options that we might want, they will go here...
	id_ = tag->getOption<std::string>( "name" );
	parse_tag( tag, data );
}

std::string const &
ConstraintGenerator::id() const
{
	return id_;
}

void
ConstraintGenerator::set_id( std::string const & id )
{
	id_ = id;
}

std::string const &
ConstraintGenerator::class_name() const
{
	return class_name_;
}

} //namespace constraint_generator
} //namespace protocols
