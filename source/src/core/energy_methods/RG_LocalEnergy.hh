// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RG_LocalEnergy.hh
/// @brief  only calculates RG the length of the repeat. Uses code from RG_Energy_fast
/// @author TJ Brunette


#ifndef INCLUDED_core_energy_methods_RG_LocalEnergy_hh
#define INCLUDED_core_energy_methods_RG_LocalEnergy_hh


// Package headers
#include <core/energy_methods/RG_Energy_Fast.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace energy_methods {



class RG_Local_MinData:  public RG_MinData {
public:
	RG_Local_MinData():RG_MinData(){}
};

typedef utility::pointer::shared_ptr< RG_Local_MinData > RG_Local_MinDataOP;

class RG_LocalEnergy: public RG_Energy_Fast  {

public:
	/// @brief Defines a center of mass based RG calculation that is O(n) rather
	/// than O(n^2).
	RG_LocalEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;


	core::Real
	calculate_rg_score( pose::Pose const & pose ) const;

	core::Real
	calculate_rg_score(
		pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues) const;

	// derivatives
	void setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sf) const override;

	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:
	RG_Local_MinData const & mindata_from_pose( pose::Pose const & ) const;
	RG_Local_MinData & nonconst_mindata_from_pose( pose::Pose & ) const;
	core::Size lastRes_;
	core::Size firstRes_;
};


}
}

#endif // INCLUDED_core_energy_methods_RG_Energy_Fast_HH
