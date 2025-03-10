// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_DataInfo.fwd.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_rna_RNA_DataInfo_fwd_hh
#define INCLUDED_core_pose_rna_RNA_DataInfo_fwd_hh

// C++

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

//Auto Headers


namespace core {
namespace pose {
namespace rna {

enum RNA_ReactivityType {
	NO_REACTIVITY,
	DMS,
	CMCT,
	SHAPE_1M7,
	HRF
};

class RNA_DataInfo;
typedef utility::pointer::shared_ptr< RNA_DataInfo > RNA_DataInfoOP;
typedef utility::pointer::weak_ptr< RNA_DataInfo > RNA_DataInfoAP;

class RNA_Reactivity;
typedef utility::pointer::shared_ptr< RNA_Reactivity > RNA_ReactivityOP;
typedef utility::pointer::weak_ptr< RNA_Reactivity > RNA_ReactivityAP;

} //rna
} //pose
} //core

#endif
