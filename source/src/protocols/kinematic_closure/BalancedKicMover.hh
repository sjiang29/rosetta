// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_BalancedKicMover_HH
#define INCLUDED_protocols_kinematic_closure_BalancedKicMover_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/BalancedKicMover.fwd.hh>
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>
#include <protocols/moves/MoverCreator.hh>

// Core headers
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Protocol headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/LoopMover.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/noncopyable.hpp>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace kinematic_closure {

/// @brief Make a kinematic closure move that obeys detailed balance.
///
/// @details This class is very conceptually similar to KicMover, so check
/// out its documentation for a general overview of the kinematic closure
/// algorithm.  Here I will just highlight some details associated with making
/// a balanced version of the move.  Detailed balance is a useful property,
/// because it allows a Monte Carlo simulation to recapitulate ensembles with
/// correct equilibrium populations (so long as sampling is good, of course).
/// There are a two reasons why the standard KicMover algorithm does not obey
/// detailed balance.  The first is that the geometry of the closure move
/// itself introduces some inherent bias which has to be explicitly canceled
/// out.  The second is that care needs to be taken to perturb the non-pivot
/// torsions in a way that also obeys detailed balance, as well.
///
/// The add_perturber() method works much like it does in KicMover.  The only
/// conceptual difference is that when the added perturbers are used internally
/// within apply(), perturbers::Perturber::perturb_with_balance() is called
/// instead of perturbers::Perturber::perturb().  This makes it easy to make
/// variants of the perturber algorithms which obey detailed balance.  The
/// set_pivot_picker() method is no different from the KicMover version.

class BalancedKicMover : public protocols::canonical_sampling::ThermodynamicMover {

public:

	/// @brief Default constructor.
	BalancedKicMover();

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Default destructor.
	~BalancedKicMover() override;

public:

	/// @copydoc KicMover::apply
	void apply(Pose & pose) override;

	/// @copydoc KicMover::get_name
	std::string get_name() const override;

	static std::string mover_name();

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_balancedKIC_mover( utility::tag::XMLSchemaDefinition & xsd );

public:

	/// @copydoc KicMover::set_loop
	void set_loop(Loop const & loop);

	/// @copydoc KicMover::set_loops
	void set_loops(protocols::loops::LoopsOP const loops);

	/// @copydoc KicMover::add_perturber
	void add_perturber(perturbers::PerturberOP perturber);

	/// @copydoc KicMover::set_pivot_picker
	void set_pivot_picker(pivot_pickers::PivotPickerOP picker);

	/// @brief Set the residue_selector to use
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Get the residue selector that this mover uses
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const{ return residue_selector_; }

public:

	/// @brief Return true, because this mover always obeys detailed balance.
	bool preserve_detailed_balance() const override { return true; }

	/// @brief This mover always obeys detailed balance, so this is a no-op.
	void set_preserve_detailed_balance(bool) override {}

	/// @details Right now the proposal probabilities are balanced internally, so
	/// this ratio will always be unity.  This could change eventually, though.
	Real last_proposal_density_ratio() override { return 1; }

	/// @brief Indicate that each torsion in the loop may take on any value.
	utility::vector1<core::id::TorsionID_Range> torsion_id_ranges(Pose & pose) override;

public:

	/// @brief Pick a solution in a way that cancels out the geometrical bias of
	/// the kinematic closure algorithm.
	static ClosureSolutionCOP pick_solution(
		SolutionList const & unperturbed_solutions,
		SolutionList const & perturbed_solutions);


	/// @brief Return true if the given solution is the same as the input pose.
	/// This allows for a more accurate reporting of Monte Carlo statistics.
	static bool is_solution_trivial(
		ClosureProblemCOP problem,
		ClosureSolutionCOP solution,
		SolutionList const & unperturbed_solutions,
		SolutionList const & perturbed_solutions);

private:
	/// Methods
	void init_from_options();

private:
	bool is_fold_tree_stale_;
	bool selector_on_;
	std::string loops_file_;
	protocols::loops::Loop loop_;
	protocols::loops::LoopsOP loops_;
	perturbers::PerturberSetOP perturbers_;
	pivot_pickers::PivotPickerOP pivot_picker_;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

};

}
}

#endif
