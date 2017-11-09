/*
 *  histentry_cs.cpp
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

/**
 * \file histentry_cs.cpp
 * Implementation of archiving_node to record and manage spike history
 * \author Jesus Garrido
 * \date october 2017
 * \note moved to separate file to avoid circular inclusion in node.h
 */

#include "histentry_cos.h"

// member functions of histentry

mynest::histentry_cos::histentry_cos( double t, double cos2, double sin2, double cossin, size_t access_counter )
  : t_( t )
  , cos2_( cos2 )
  , sin2_( sin2 )
  , cossin_( cossin )
  , access_counter_( access_counter )
{
}
