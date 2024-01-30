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
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

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


    return 0;
}