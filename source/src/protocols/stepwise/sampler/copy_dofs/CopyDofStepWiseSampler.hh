// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/copy_dofs/CopyDofStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_copy_dofs_CopyDofStepWiseSampler_HH
#define INCLUDED_protocols_sampler_copy_dofs_CopyDofStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/copy_dofs/CopyDofStepWiseSampler.fwd.hh>
#include <protocols/simple_moves/CopyDofMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh> // AUTO IWYU For vector1
#include <map> // AUTO IWYU For map

#ifdef WIN32
#include <protocols/simple_moves/CopyDofMover.hh>
#endif

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

class CopyDofStepWiseSampler: public StepWiseSamplerSized {

public:

	//constructor
	CopyDofStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
		std::map< core::Size, core::Size > const & res_map,
		core::pose::Pose const & starting_pose );

	//constructor
	CopyDofStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
		std::map< core::Size, core::Size > const & res_map );


public:

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const override{ return copy_dof_movers_.size(); }

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose&, core::Size const ) override;

	/// @brief Name of the class
	std::string get_name() const override { return "CopyDofStepWiseSampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::COPY_DOF; }


protected:
	utility::vector1< simple_moves::CopyDofMoverOP > copy_dof_movers_;
	utility::vector1< core::pose::PoseOP > pose_list_;


};

} //copy_dofs
} //sampler
} //stepwise
} //protocols

#endif
