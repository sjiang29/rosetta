#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>

int main(int argc, char ** argv) {

    std::cout << "Hello World!" << std::endl;
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
        std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }

    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
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


    protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo( *mypose, score, 25);

    for(int i = 0; i < 5; i++){
        double uniform_random_number = numeric::random::uniform();
        core::Size randres = uniform_random_number * (mypose->total_residue() + 1);//… code here to pick the index of a random residue in the Pose
        core::Real pert1 = numeric::random::gaussian();//… code here to get a random number
        core::Real pert2 = numeric::random::gaussian();//… code here to get another random number
        core::Real orig_phi = mypose->phi( randres );
        core::Real orig_psi = mypose->psi( randres );
        mypose->set_phi( randres, orig_phi + pert1 );
        mypose->set_psi( randres, orig_psi + pert2 );

        mc.boltzmann( *mypose, score);

    }



    return 0;
}