// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/TreeBuilder.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_TREEBUILDER_FWD_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_TREEBUILDER_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace nonlocal {

class TreeBuilder;
typedef utility::pointer::shared_ptr<TreeBuilder> TreeBuilderOP;
typedef utility::pointer::shared_ptr<TreeBuilder const> TreeBuilderCOP;

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_TREEBUILDER_FWD_HH_
