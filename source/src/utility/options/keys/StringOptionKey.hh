// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/StringOptionKey.hh
/// @brief  Automatic hidden index key for string options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_StringOptionKey_hh
#define INCLUDED_utility_options_keys_StringOptionKey_hh


// Unit headers
#include <utility/options/keys/StringOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/StringOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace utility {
namespace options {


/// @brief Automatic hidden index key for string options
class StringOptionKey :
	public ScalarOptionKey
{


private: // Types


	typedef  ScalarOptionKey  Super;


private: // Friends


	friend class utility::keys::SmallKeyVector< StringOptionKey, StringOption >;


public: // Creation


	/// @brief Default constructor
	inline
	StringOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	StringOptionKey(
		StringOptionKey const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Key constructor
	inline
	explicit
	StringOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	StringOptionKey(
		Key const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Identifier constructor
	inline
	explicit
	StringOptionKey(
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Clone this
	StringOptionKey *
	clone() const override
	{
		return new StringOptionKey( *this );
	}


	/// @brief Destructor
	inline
	~StringOptionKey() override {}


public: // Assignment


	/// @brief Key assignment
	inline
	StringOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // StringOptionKey


} // namespace options
} // namespace utility


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( utility_options_keys_StringOptionKey )
#endif // SERIALIZATION


#endif // INCLUDED_utility_options_keys_StringOptionKey_HH
