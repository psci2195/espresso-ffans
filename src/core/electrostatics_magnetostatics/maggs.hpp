/*
 Copyright (C) 2010,2011 Florian Fahrenberger
 Copyright (C) 2010-2018 The ESPResSo project
 Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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

/** \file
 *  Maxwell Equations Molecular Dynamics (MEMD) method for electrostatic
 *  interactions.
 *
 *  We use a local update scheme to propagate artificial B-fields on a
 *  lattice in the system. In principal, the algorithm simulates full
 *  electrodynamics, but with a tunable speed of light.
 *
 *  The method is very usable for large particle numbers or highly
 *  parallel architectures, since it is local and scales linearly.
 *  It is not suited for high-precision calculation of forces, since
 *  the simple interpolation scheme produces errors in the order of
 *  10^-5 in the force.
 *
 *  The chosen mesh should roughly be of the size of the particles.
 *
 *  Further reading on the algorithm:
 *  <ul>
 *  <li> I. Pasichnyk and B. Dunweg, Coulomb interaction via local dynamics: a
 * molecular-dynamics algorithm. J. Phys: Condens. Matter, 16 ,p. 1399-4020,
 * (2004).
 *  </ul>
 *
 */

#ifndef _MAGGS_H
#define _MAGGS_H

#include "config.hpp"

#ifdef ELECTROSTATICS

/** \name External structure */
/*@{*/

/** Global system information for MEMD algorithm.
*/
typedef struct {
  int finite_epsilon_flag;
  int adaptive_flag;
  double scaling;
  double epsilon_infty;
  /** speed of light parameter */
  double f_mass; // = 1/c^2
  /** inverse of square root of f_mass */
  double invsqrt_f_mass;
  double prefactor;
  /** prefactor to convert field to force */
  double pref1;
  /** mesh size in one dimension */
  int mesh;
  /* = 1/a = mesh / box_length */
  double inva;
  /** size of mesh cube */
  double a;
  /** self energy coefficients of the system */
  double alpha[8][8];
} MAGGS_struct;
extern MAGGS_struct maggs;
/*@}*/

/*****************************/
/** \name External functions */
/*****************************/

/*@{*/

/** Initialization function, parse command and set parameters.
    Called from initialize.cpp
*/
void maggs_init();

/** Set the main parameters for the algorithm.
    @param prefactor   Electrostatics prefactor for the system
    @param f_mass    parameter to tune the speed of light (1/c^2)
    @param mesh      Mesh size in one dimension
    @param finite_epsilon_flag whether to do epsilon-at-infinity-correction
    @param epsilon_infty epsilon-at-infinity
 */
int maggs_set_parameters(double prefactor, double f_mass, int mesh,
                         int finite_epsilon_flag, double epsilon_infty);

/** Get lattice size in one dimension.
 @return mesh in 1D
 */
int maggs_get_mesh_1D();

/** Set permittivity for single lattice links.
 @param node_x              index of the node in x direction
 @param node_y              index of the node in y direction
 @param node_z              index of the node in z direction
 @param direction           direction in which the link points from the node. 0
 is for x, 1 is for y, 2 is for z
 @param relative_epsilon    permittivity to set, relative to the background
 permittivity set by the electrostatics prefactor
 */
double maggs_set_permittivity(int node_x, int node_y, int node_z, int direction,
                              double relative_epsilon);

/** Set adaptive permittivity flag.
 @param scaling             scaling of the volumetric formula for salt dependent
 permittivity
 */
int maggs_set_adaptive_flag(double scaling);

/** Propagate the B-field in the system.
    Called TWICE from \ref integrate.cpp with timestep dt/2 to ensure
   time-reversibility of the integrator.
    @param dt Timestep for which to propagate the field.
*/
void maggs_propagate_B_field(double dt);

/** Calculate the forces on all particles. Writes the result directly into the
 * force pointer. */
void maggs_calc_forces();

/** Get the electric energy of the system as return value */
double maggs_electric_energy();
/** Get the magnetic energy of the artificial transversal field component as
 * return value */
double maggs_magnetic_energy();

/** Count the number of charges in the whole system. */
int maggs_count_charged_particles();

/** Clean up, free memory. Not used at the moment. */
void maggs_exit();

/*@}*/

#endif
#endif
