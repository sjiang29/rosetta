// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddPDBInfoMover.hh
/// @author Jonathan Weinstein jonathan.weinstein@weizmann.ac.il

#ifndef INCLUDED_protocols_simple_moves_AddPDBInfoMover_hh
#define INCLUDED_protocols_simple_moves_AddPDBInfoMover_hh

#include <protocols/simple_moves/AddPDBInfoMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>



namespace protocols {
namespace simple_moves {

class AddPDBInfoMover : public protocols::moves::Mover
{
public:
	AddPDBInfoMover();
	~AddPDBInfoMover() override;
	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
};


} // simple_moves
} // protocols


#endif /*INCLUDED_protocols_moves_AddPDBInfoMover_HH*/
