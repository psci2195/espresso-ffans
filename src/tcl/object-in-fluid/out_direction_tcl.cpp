/*
  Copyright (C) 2012,2013 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "tcl/parser.hpp"
#include "object-in-fluid/out_direction.hpp"
#include "out_direction_tcl.hpp"
#include "interaction_data.hpp"

#ifdef MEMBRANE_COLLISION
// parse parameters for the out_direction
int tclcommand_inter_parse_oif_out_direction(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  
  CHECK_VALUE(out_direction_set_params(bond_type), "bond type must be nonnegative");
}

int tclprint_to_result_oifoutdirectionIA(Tcl_Interp *interp, Bonded_ia_parameters *params){

// no parameters to be added here
    return (TCL_OK);
}
#endif

