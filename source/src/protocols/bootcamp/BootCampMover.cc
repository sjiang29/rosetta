// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief a subclass of mover
/// @author sjiang29 (69825128+sjiang29@users.noreply.github.com)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( "BootCampMover" )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& mypose){
            core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
            core::Real score = sfxn->score( *mypose );

            std::cout << "The score of sfxn is:" << score << std::endl;

            //double uniform_random_number = numeric::random::uniform();
            //core::Size randres = uniform_random_number * (mypose->total_residue() + 1);//… code here to pick the index of a random residue in the Pose
            //core::Real pert1 = numeric::random::gaussian();//… code here to get a random number
            //core::Real pert2 = numeric::random::gaussian();//… code here to get another random number
            //core::Real orig_phi = mypose->phi( randres );
            //core::Real orig_psi = mypose->psi( randres );
            //mypose->set_phi( randres, orig_phi + pert1 );
            //mypose->set_psi( randres, orig_psi + pert2 );


            protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo( *mypose, *sfxn, 1);
            protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
            the_observer->pymol().apply( *mypose);


            int accept_count = 0;
            int iteration_count = 200;
            for(int i = 0; i < iteration_count; i++){
                double uniform_random_number = numeric::random::uniform(); // get a random number in a range [0, 1]
                core::Size randres = uniform_random_number * mypose->total_residue() + 1;//… code here to pick the index of a random residue in the Pose
                core::Real pert1 = numeric::random::gaussian();//… code here to get a random number [-inf, inf]
                core::Real pert2 = numeric::random::gaussian();//… code here to get another random number
                core::Real orig_phi = mypose->phi( randres );
                core::Real orig_psi = mypose->psi( randres );
                mypose->set_phi( randres, orig_phi + pert1 );
                mypose->set_psi( randres, orig_psi + pert2 );

                if (mc.boltzmann( *mypose)){
                    accept_count++;
                }
            }

            double accept_rate = static_cast<double>(accept_count) / iteration_count;
            double rejection_rate = 1 - accept_rate;

            std::cout << "The accept rate is:" << accept_rate << std::endl;


            core::pack::task::PackerTaskOP repack_task =
                    core::pack::task::TaskFactory::create_packer_task( *mypose );
            repack_task->restrict_to_repacking();
            core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

            core::kinematics::MoveMap mm;
            mm.set_bb( true );
            mm.set_chi( true );

            core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
            core::optimization::AtomTreeMinimizer atm;
            atm.run( *mypose, mm, *sfxn, min_opts );

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "a subclass of mover", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns sjiang29 as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"sjiang29",
		"TODO: institution",
		"69825128+sjiang29@users.noreply.github.com",
		"Wrote the BootCampMover."
		)
	);
}

void
BootCampMover::set_num_iterations_(core::size num_iterations)  { num_iterations_ = num_iterations; }

core::size
BootCampMover::get_num_iterations_() const { return num_iterations_; }


void BootCampMover::score_function( ScoreFunctionOP sf )
{
    runtime_assert( sf != nullptr );
    sfxn_ = sf;
}

ScoreFunctionOP
BootCampMover::get_score_function_() const { return sfxn_; }




////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
