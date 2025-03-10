// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/chunk_library_inputters/chunk_library_inputter_schemas.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_chunk_library_inputters_chunk_library_inputter_schemas_HH
#define INCLUDED_protocols_jd3_chunk_library_inputters_chunk_library_inputter_schemas_HH

// Project headers

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace chunk_library_inputters {

/// @brief Define the XML schema definition for a ChunkLibraryInputter that contains
/// no subtags but may contain any number of attributes (aka options).
void
chunk_library_inputter_xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & inputter_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);


}  // namespace chunk_library_inputters
}  // namespace jd3
}  // namespace protocols

#endif
