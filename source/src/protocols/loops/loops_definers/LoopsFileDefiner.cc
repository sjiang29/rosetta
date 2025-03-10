// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsFileDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsFileDefiner.hh>

// Package headers
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_definers/util.hh>

// Project Headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>


using std::string;
using std::endl;
using std::stringstream;
using utility::tag::TagCOP;
using basic::datacache::DataMap;
using core::pose::Pose;


namespace protocols {
namespace loops {
namespace loops_definers {

LoopsFileDefiner::LoopsFileDefiner() = default;

LoopsFileDefiner::~LoopsFileDefiner() = default;

LoopsFileDefiner::LoopsFileDefiner(LoopsFileDefiner const & /*src*/) = default;

/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsFileDefiner::clone(
) const {
	return utility::pointer::make_shared< LoopsFileDefiner >(*this);
}

/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsFileDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &
) {

	if ( !tag->hasOption("name") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,
			"Unable to create unnamed LoopsDefiner (type: LoopsFile)" );
	}
	string const loops_name(tag->getOption<string>("name"));


	string filename;
	if ( tag->hasOption("filename") ) {
		filename = tag->getOption<string>("filename");
	} else {
		stringstream err_msg;
		err_msg << "Tag with name '" << loops_name << "' does not have the expected 'filename' field." << endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, err_msg.str());
	}

	LoopsFileIO loops_file_io;
	lfd_ = loops_file_io.read_loop_file( filename );
}

SerializedLoopList
LoopsFileDefiner::apply(
	Pose const & pose
) {
	runtime_assert( lfd_ != nullptr );
	return lfd_->resolve_as_serialized_loops( pose );
}

std::string LoopsFileDefiner::class_name()
{
	return "LoopsFile";
}

void LoopsFileDefiner::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using Attr = XMLSchemaAttribute;

	AttributeList attributes;
	attributes + Attr::required_attribute( "filename", xs_string, "The file from which the loops should be read" );

	xsd_type_definition_w_attributes( xsd, class_name(), "Define a set of loops reading them in from the provided file", attributes );
}

} //namespace
} //namespace
} //namespace
