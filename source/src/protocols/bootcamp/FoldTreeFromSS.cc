//
// Created by Shan Jiang on 1/31/24.
//

#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>
#include <protocols/loops/loop>


namespace protocols{
namespace  bootcamp{


    FoldTreeFromSS::FoldTreeFromSS(std::string const & ssstring){
        int ss_start = 1;
        int ss_end = ss.size();

        utility::vector1< std::pair< core::Size, core::Size > > pairs = identify_secondary_structure_spans(ss);

        core::kinematics::FoldTree ft;

        std::pair< core::Size, core::Size > first_ele = pairs[1];
        core::Size start1 = first_ele.first;
        core::Size ending1 = first_ele.second;

        core::Size mid1 = get_middle(start1, ending1);

        int jump_count = 1;
        for( unsigned long i = 2; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];
            std::pair< core::Size, core::Size > prev_ele = pairs[i-1];

            core::Size start2 = ele.first;
            core::Size ending2 = ele.second;
            core::Size mid2 = get_middle(start2, ending2);
            // jump from mid residue of first ele to mid residue of all other elements
            ft.add_edge( mid1,mid2,jump_count++);

            core::Size ending3 = prev_ele.second;
            // jump from mid residue to mid residue of every inter-secondary-gap
            core::Size mid3 = get_middle(start2, ending3);
            ft.add_edge( mid1, mid3, jump_count++);
        }

        for( unsigned long i = 1; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];

            core::Size start2 = ele.first;
            core::Size ending2 = ele.second;
            core::Size mid2 = get_middle(start2, ending2);

            if(i == 1){
                ft.add_edge( mid2,ending2, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,ss_start, core::kinematics::Edge::PEPTIDE);
            } else if (i == pairs.size()) {
                ft.add_edge( mid2,ss_end, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,start2, core::kinematics::Edge::PEPTIDE);
            } else {
                ft.add_edge( mid2,ending2, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,start2, core::kinematics::Edge::PEPTIDE);
            }


        }

        for( unsigned long i = 2; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];
            std::pair< core::Size, core::Size > prev_ele = pairs[i-1];

            core::Size curr_start = ele.first;
            //core::Size ending2 = ele.second;
            //core::Size mid2 = get_middle(start2, ending2);
            // jump from mid residue of first ele to mid residue of all other elements
            //ft.add_edge( mid1,mid2,jump_count++);

            core::Size prev_ending = prev_ele.second;
            // jump from mid residue to mid residue of every inter-secondary-gap
            core::Size mid = get_middle(curr_start, prev_ending);
            ft.add_edge( mid, prev_ending, core::kinematics::Edge::PEPTIDE);
            ft.add_edge( mid, curr_start, core::kinematics::Edge::PEPTIDE);

        }

    }

    FoldTreeFromSS::~FoldTreeFromSS() = default;

    core::kinematics::FoldTree fold_tree_from_dssp_string(std::string ss){

        int ss_start = 1;
        int ss_end = ss.size();

        utility::vector1< std::pair< core::Size, core::Size > > pairs = identify_secondary_structure_spans(ss);

        core::kinematics::FoldTree ft;

        std::pair< core::Size, core::Size > first_ele = pairs[1];
        core::Size start1 = first_ele.first;
        core::Size ending1 = first_ele.second;

        core::Size mid1 = get_middle(start1, ending1);

        int jump_count = 1;
        for( unsigned long i = 2; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];
            std::pair< core::Size, core::Size > prev_ele = pairs[i-1];

            core::Size start2 = ele.first;
            core::Size ending2 = ele.second;
            core::Size mid2 = get_middle(start2, ending2);
            // jump from mid residue of first ele to mid residue of all other elements
            ft.add_edge( mid1,mid2,jump_count++);

            core::Size ending3 = prev_ele.second;
            // jump from mid residue to mid residue of every inter-secondary-gap
            core::Size mid3 = get_middle(start2, ending3);
            ft.add_edge( mid1, mid3, jump_count++);
        }

        for( unsigned long i = 1; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];

            core::Size start2 = ele.first;
            core::Size ending2 = ele.second;
            core::Size mid2 = get_middle(start2, ending2);

            if(i == 1){
                ft.add_edge( mid2,ending2, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,ss_start, core::kinematics::Edge::PEPTIDE);
            } else if (i == pairs.size()) {
                ft.add_edge( mid2,ss_end, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,start2, core::kinematics::Edge::PEPTIDE);
            } else {
                ft.add_edge( mid2,ending2, core::kinematics::Edge::PEPTIDE);
                ft.add_edge( mid2,start2, core::kinematics::Edge::PEPTIDE);
            }


        }

        for( unsigned long i = 2; i <= pairs.size(); i++) {
            std::pair< core::Size, core::Size > ele = pairs[i];
            std::pair< core::Size, core::Size > prev_ele = pairs[i-1];

            core::Size curr_start = ele.first;
            //core::Size ending2 = ele.second;
            //core::Size mid2 = get_middle(start2, ending2);
            // jump from mid residue of first ele to mid residue of all other elements
            //ft.add_edge( mid1,mid2,jump_count++);

            core::Size prev_ending = prev_ele.second;
            // jump from mid residue to mid residue of every inter-secondary-gap
            core::Size mid = get_middle(curr_start, prev_ending);
            ft.add_edge( mid, prev_ending, core::kinematics::Edge::PEPTIDE);
            ft.add_edge( mid, curr_start, core::kinematics::Edge::PEPTIDE);

        }
        return ft;
    }
}
}