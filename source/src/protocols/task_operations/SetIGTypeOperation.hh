// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/SetIGTypeOperation.hh
///
/// @brief Task operation to set interaction graph type (linear memory, lazy or double lazy)
/// @author Sagar Khare


#ifndef INCLUDED_protocols_task_operations_SetIGTypeOperation_HH
#define INCLUDED_protocols_task_operations_SetIGTypeOperation_HH

// Unit Headers
#include <protocols/task_operations/SetIGTypeOperation.fwd.hh>

// Core Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace task_operations {

class SetIGTypeOperation : public core::pack::task::operation::TaskOperation
{
public:
	SetIGTypeOperation();

	SetIGTypeOperation( bool lin_mem, bool lazy, bool double_lazy, bool precompute );

	~SetIGTypeOperation() override;

	core::pack::task::operation::TaskOperationOP clone() const override;

	void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const override;

	void parse_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SetIGType"; }

private:

	bool lin_mem_,lazy_,double_lazy_,precompute_;

};

} //task_operations
} //protocols

#endif
