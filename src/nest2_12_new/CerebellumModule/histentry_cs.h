/*
 *  histentry_cs.h
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
 * \file histentry_cs.h
 * Part of definition of Archiving_Node which is capable of
 * recording and managing a spike history.
 * \author Jesus Garrido
 * \date october 2017
 */

#ifndef HISTENTRY_CS_H
#define HISTENTRY_CS_H

// Includes from nestkernel:
#include "nest_types.h"

namespace mynest
{

// entry in the spiking history
class histentry_cs
{
public:
  histentry_cs( double t, size_t access_counter );

  double t_;              //!< point in time when spike occurred (in ms)
  //! how often this entry was accessed (to enable removal, once read by all
  //! neurons which need it)
  size_t access_counter_;
};
}

#endif
