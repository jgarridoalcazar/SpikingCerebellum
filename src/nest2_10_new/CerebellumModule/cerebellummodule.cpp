/*
 *  cerebellummodule.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"
#include "connector_model_impl.h"
#include "target_identifier.h"

// include headers with your own stuff
#include "cerebellummodule.h"
#include "stdp_sin_connection.h"
#include "iaf_cond_exp_cs.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module.
 * You make the module known to the loader by defining an instance of your
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::CerebellumModule cerebellummodule_LTX_mod;

// -- DynModule functions ------------------------------------------------------

mynest::CerebellumModule::CerebellumModule()
{
#ifdef LINKED_MODULE
  // register this module at the dynamic loader
  // this is needed to allow for linking in this module at compile time
  // all registered modules will be initialized by the main app's dynamic loader
  nest::DynamicLoaderModule::registerLinkedModule( this );
#endif
}

mynest::CerebellumModule::~CerebellumModule()
{
}

const std::string
mynest::CerebellumModule::name( void ) const
{
  return std::string( "Cerebellum Module" ); // Return name of the module
}

const std::string
mynest::CerebellumModule::commandstring( void ) const
{
  // Instruct the interpreter to load mymodule-init.sli
  return std::string( "(cerebellummodule-init) run" );
}

void
mynest::CerebellumModule::init( SLIInterpreter* i )
{
  /* Register a neuron or device model.
     Give node type as template argument and the name as second argument.
     The first argument is always a reference to the network.
  */
  nest::register_model< mynest::iaf_cond_exp_cs >( nest::NestModule::get_network(), "iaf_cond_exp_cs" );

  /* Register a synapse type.
     Give synapse type as template argument and the name as second argument.
     The first argument is always a reference to the network.

     There are two choices for the template argument:
         - nest::TargetIdentifierPtrRport
         - nest::TargetIdentifierIndex
     The first is the standard and you should usually stick to it.
     nest::TargetIdentifierIndex reduces the memory requirement of synapses
     even further, but limits the number of available rports. Please see
     Kunkel et al, Front Neurofinfom 8:78 (2014), Sec 3.3.2, for details.
  */
  nest::register_connection_model< mynest::STDPSinConnection< nest::TargetIdentifierPtrRport > >(
    nest::NestModule::get_network(), "stdp_sin_synapse" );

} // MyModule::init()
