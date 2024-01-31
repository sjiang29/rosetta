//
// Created by Shan Jiang on 2/1/24.
//


#ifndef INCLUDED_protocols_bootcamp_fold_tree_from_ss_HH
#define INCLUDED_protocols_bootcamp_fold_tree_from_ss_HH

#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace bootcamp  {
    utility::vector1< std::pair< core::Size, core::Size > >
    identify_secondary_structure_spans( std::string const & ss_string );

    core::size get_middle(int start, int end);

    core::kinematics::FoldTree fold_tree_from_dssp_string(std::string ss);

    core::kinematics::FoldTree fold_tree_from_dssp_string(core::pose::Pose const & init_pose);

}

}

#endif //INCLUDED_protocols_bootcamp_fold_tree_from_ss_HH