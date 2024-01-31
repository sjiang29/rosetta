//
// Created by Shan Jiang on 1/31/24.
//

#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>

utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
    utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
    core::Size strand_start = -1;
    for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
        if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
            if ( int( strand_start ) == -1 ) {
                strand_start = ii;
            } else if ( ss_string[ii] != ss_string[strand_start] ) {
                ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                strand_start = ii;
            }
        } else {
            if ( int( strand_start ) != -1 ) {
                ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
                strand_start = -1;
            }
        }
    }
    if ( int( strand_start ) != -1 ) {
        // last residue was part of a ss-eleemnt
        ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
    }
    for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
        std::cout << "SS Element " << ii << " from residue "
                  << ss_boundaries[ ii ].first << " to "
                  << ss_boundaries[ ii ].second << std::endl;
    }
    return ss_boundaries;
}

core::size get_middle(int start, int end) {
    return (start + end) / 2;
}

core::kinematics::FoldTree fold_tree_from_dssp_string(std::string ss){

    utility::vector1< std::pair< core::Size, core::Size > > pairs = identify_secondary_structure_spans(ss);
    core::kinematics::FoldTree ft;

    std::pair< core::Size, core::Size > first_ele = pairs[1];
    core::Size start1 = first_ele.first;
    core::Size ending1 = first_ele.second;

    core::Size mid1 = get_middle(start1, ending1);

    int jump_count = 1;
    for( unsigned long i = 2; i != pairs.size(); i++) {
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

    for( unsigned long i = 1; i != pairs.size(); i++) {
        std::pair< core::Size, core::Size > ele = pairs[i];

        core::Size start2 = ele.first;
        core::Size ending2 = ele.second;
        core::Size mid2 = get_middle(start2, ending2);

        ft.add_edge( mid2,ending2, core::kinematics::Edge::PEPTIDE);
        ft.add_edge( mid2,start2, core::kinematics::Edge::PEPTIDE);

    }
    return ft;
}
