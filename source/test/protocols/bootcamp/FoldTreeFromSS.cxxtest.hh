// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>



#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers

/// Project headers
#include <core/types.hh>

// C++ headers
#include <core/kinematics/FoldTree.hh>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>





// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	void test_hello_world() {
        TS_ASSERT(true);
    }

    void test_identify_secondary_structure_spans() {
        std::string s = "   EEEEE   ";
        utility::vector1< std::pair< core::Size, core::Size > > res = identify_secondary_structure_spans(s);
        // everything is 1-based
        std::pair< core::Size, core::Size > p = res[1];
        core::Size start = p.first;
        core::Size ending = p.second;

        TS_ASSERT_EQUALS(4, start);
    }

    void find_ones_block( utility::vector1< int > const & bitstring ) {
        int start = 0;
        for ( uint ii = 1; ii <= bitstring.size(); ++ii ) {
            // std::cout << "ii: " << ii << " " << bitstring[ ii ] << std::endl;
            if ( start != 0 ) {
                if ( bitstring[ ii ] != 1 ) {
                    std::cout << "start: " << start << " stop: " << ii-1 << std::endl;
                    start = 0;
                }
            } else {
                if ( bitstring[ ii ] == 1 ) {
                    start = ii;
                }
            }
        }
    }

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

    core::Size get_middle(int start, int end) {
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






};
