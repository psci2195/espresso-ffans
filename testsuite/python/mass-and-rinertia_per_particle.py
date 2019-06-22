# Copyright (C) 2010-2018 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import uniform
import espressomd
import math
import random


@ut.skipIf(not espressomd.has_features(["MASS",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "LANGEVIN_PER_PARTICLE"]),
           "Features not available, skipping test!")
class ThermoTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))
    gamma_global_rot = np.zeros((3))
    # Test ranges    
    gamma_min = 5.
    gamma_max = 10.

    # Particle properties
    mass = 0.0
    J = 0.0, 0.0, 0.0

    # Per-particle type parameters.
    # 2 different langevin parameters for particles.
    kT_p = np.zeros((2))
    # gamma_tran/gamma_rot matrix: [2 kinds of particles] x [3 dimensions X Y Z]
    # These matrices are assigning per-particle in corresponding test cases.
    gamma_tran_p = np.zeros((2, 3))
    gamma_rot_p = np.zeros((2, 3))

    # These variables will take the values to compare with.
    # Depending on the test case following matrices either equals to the previous
    # or the global corresponding parameters. The corresponding setting effect is an essence of
    # all the test cases' differentiation here.
    halfkT_p_validate = np.zeros((2))
    gamma_tran_p_validate = np.zeros((2, 3))
    gamma_rot_p_validate = np.zeros((2, 3))
    # Diffusivity
    D_tran_p_validate = np.zeros((2, 3))

    @classmethod
    def setUpClass(cls):
        # Handle a random generator seeding
        seed1 = 15
        np.random.seed(seed1)
        seed2 = 42
        # The Espresso system configuration
        cls.system.seed = [s * seed2 for s in range(cls.system.cell_system.get_state()["n_nodes"])]
        cls.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        cls.system.cell_system.skin = 5.0

    def setUp(self):
        self.system.time = 0.0
        self.system.part.clear()
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.system.thermostat.turn_off()
    
    def set_initial_cond(self):
        """
        Set all the particles to zero coordinates and velocities; same for time.
        The quaternion is set to default value.

        """
        system = self.system
        system.time = 0.0
        system.part[:].pos = np.zeros((3))
        system.part[:].v = np.zeros((3))
        system.part[:].omega_body = np.zeros((3))
        system.part[:].quat = np.array([1., 0., 0., 0.])

    def set_langevin_global_defaults(self):
        """
        Setup the expected NVT thermostat viscous friction parameters.

        """

        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k,:] = self.gamma_global[:]
            self.gamma_rot_p_validate[k,:] = self.gamma_global[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def set_langevin_global_defaults_rot_differ(self):
        """
        Setup the expected NVT thermostat viscous friction parameters
        with a rotation-specific gamma.

        """
        # Global NVT thermostat parameters are assigning by default
        for k in range(2):
            self.gamma_tran_p_validate[k,:] = self.gamma_global[:]
            self.gamma_rot_p_validate[k,:] = self.gamma_global_rot[:]
            self.halfkT_p_validate[k] = self.kT / 2.0

    def dissipation_param_setup(self, n):
        """
        Setup the parameters for the following dissipation
        test.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        system = self.system
        # Time
        self.system.time_step = 0.007

        # Space
        box = 1.0
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        self.kT = 0.0
        # The translational gamma isotropy is required here.
        # Global gamma for tests without particle-specific gammas.
        #
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the gamma less than
        # some maximal value according to the value gamma_max.
        # Also, it cannot be very small (gamma_min), otherwise the thermalization will require
        # too much of the CPU time. Same: for all such gamma assignments throughout the test.
        #
        gamma_min = self.gamma_min
        gamma_max = self.gamma_max
        gamma_rnd = uniform(gamma_min, gamma_max)
        self.gamma_global = gamma_rnd, gamma_rnd, gamma_rnd
        # Additional test case for the specific global rotational gamma set.
        self.gamma_global_rot = uniform(gamma_min, gamma_max, 3)
        # Per-paricle values:
        self.kT_p = 0.0, 0.0
        # Either translational friction isotropy is required
        # or both translational and rotational ones.
        # Otherwise these types of motion will interfere.
        # ..Let's test both cases depending on the particle index.
        self.gamma_tran_p[0, 0] = uniform(gamma_min, gamma_max)
        self.gamma_tran_p[0, 1] = self.gamma_tran_p[0, 0]
        self.gamma_tran_p[0, 2] = self.gamma_tran_p[0, 0]
        self.gamma_rot_p[0,:] = uniform(gamma_min, gamma_max, 3)
        self.gamma_tran_p[1, 0] = uniform(gamma_min, gamma_max)
        self.gamma_tran_p[1, 1] = self.gamma_tran_p[1, 0]
        self.gamma_tran_p[1, 2] = self.gamma_tran_p[1, 0]
        self.gamma_rot_p[1, 0] = uniform(gamma_min, gamma_max)
        self.gamma_rot_p[1, 1] = self.gamma_rot_p[1, 0]
        self.gamma_rot_p[1, 2] = self.gamma_rot_p[1, 0]

        # Particles
        self.mass = 12.74
        self.J = 10.0, 10.0, 10.0

        for i in range(n):
            for k in range(2):
                ind = i + k * n
                self.system.part.add(
                    rotation=(1, 1, 1), pos=(0.0, 0.0, 0.0), id=ind)
                self.system.part[ind].v = 1.0, 1.0, 1.0
                if "ROTATION" in espressomd.features():
                    self.system.part[ind].omega_body = 1.0, 1.0, 1.0
                self.system.part[ind].mass = self.mass
                self.system.part[ind].rinertia = self.J
    
    def dissipation_viscous_drag_setup_bd(self):
        """
        Setup the specific parameters for the following dissipation
        test of the viscous drag terminal velocity stationarity.
        It is used by the BD test cases only for the moment.

        """

        system = self.system
        ## Time
        # Large time_step is OK for the BD by its definition & its benefits
        self.system.time_step = 17.0

        # Space
        # A large box is reuired due to huge steps over the space:
        # large time_step is natural for the BD.
        box = 1.0E4
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        self.kT = 0.0
        # The translational gamma isotropy is required here.
        # Global gamma for tests without particle-specific gammas.
        #
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the gamma less than
        # some maximal value according to the value max_gamma_param.
        # Also, it cannot be very small (min_gamma_param), otherwise the thermalization will require
        # too much of the CPU time. Same: for all such gamma assignments throughout the test.
        #
        gamma_rnd = uniform(self.gamma_min, self.gamma_max)
        self.gamma_global = gamma_rnd, gamma_rnd, gamma_rnd
        # Additional test case for the specific global rotational gamma set.
        gamma_rnd = uniform(self.gamma_min, self.gamma_max)
        self.gamma_global_rot = gamma_rnd, gamma_rnd, gamma_rnd
        # Per-paricle values:
        self.kT_p = 0.0, 0.0
        # Isotropy is required here for the drag tests
        for k in range(2):
            gamma_rnd_tran = uniform(self.gamma_min, self.gamma_max)
            gamma_rnd_rot = uniform(self.gamma_min, self.gamma_max)
            for j in range(3):
                self.gamma_tran_p[k, j] = gamma_rnd_tran
                self.gamma_rot_p[k, j] = gamma_rnd_rot

        # Particles
        self.mass = 12.74
        self.J = 10.0, 10.0, 10.0
        for i in range(2):
            system.part.add(rotation=(1, 1, 1), pos=(0.0, 0.0, 0.0), id=i)
            system.part[i].v = 0.0, 0.0, 0.0
            if "ROTATION" in espressomd.features():
                system.part[i].omega_body = 0.0, 0.0, 0.0
            system.part[i].mass = self.mass
            system.part[i].rinertia = self.J

    def fluctuation_dissipation_param_setup(self, n):
        """
        Setup the parameters for the following fluctuation-dissipation
        test.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """
        # Time
        self.system.time_step = 0.03

        # Space
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            # for large steps and multi-core run stability
            box = 1E2
        else:
            box = 10.0
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = uniform(1.5, 5.)
        # See the above comment regarding the gamma assignments.
        # Note: here & hereinafter specific variations in these ranges are related to
        # the test execution duration to achieve the required statistical
        # averages faster.
        gamma_min = self.gamma_min
        gamma_max = self.gamma_max
        self.gamma_global = uniform(gamma_min, gamma_max, 3)
        self.gamma_global_rot = uniform(gamma_min, gamma_max, 3)
        # Per-particle parameters
        self.kT_p = 2.5, 2.0
        for k in range(2):
            self.gamma_tran_p[k,:] = uniform(gamma_min, gamma_max, 3)
            self.gamma_rot_p[k,:] = uniform(gamma_min, gamma_max, 3)

        # Particles
        # As far as the problem characteristic time is t0 ~ mass / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the mass higher than
        # some minimal value according to the value min_mass_param.
        # Also, it is expected to test the large enough mass (max_mass_param).
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time.
        min_mass_param = 3.
        max_mass_param = 10.0
        self.mass = uniform(min_mass_param, max_mass_param)
        self.J = uniform(min_mass_param, max_mass_param, 3)
        for i in range(n):
            for k in range(2):
                ind = i + k * n
                part_pos = np.random.random(3) * box
                part_v = 0.0, 0.0, 0.0
                part_omega_body = 0.0, 0.0, 0.0
                self.system.part.add(
                    rotation=(
                        1,
                        1,
                        1),
                    id=ind,
                    mass=self.mass,
                    rinertia=self.J,
                    pos=part_pos,
                    v=part_v)
                if "ROTATION" in espressomd.features():
                    self.system.part[ind].omega_body = part_omega_body

    def check_dissipation(self, n):
        """
        Check the dissipation relations: the simple viscous decelleration test.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """
        tol = 1.3E-3
        for step in range(100):
            self.system.integrator.run(2)
            for i in range(n):
                for k in range(2):
                    ind = i + k * n
                    for j in range(3):
                        # Note: velocity is defined in the lab frame of reference
                        # while gamma_tr is defined in the body one.
                        # Hence, only isotropic gamma_tran_p_validate could be
                        # tested here.
                        self.assertLess(abs(
                            self.system.part[ind].v[j] - math.exp(- self.gamma_tran_p_validate[k, j] * self.system.time / self.mass)), tol)
                        if "ROTATION" in espressomd.features():
                            self.assertLess(abs(
                                self.system.part[ind].omega_body[j] - math.exp(- self.gamma_rot_p_validate[k, j] * self.system.time / self.J[j])), tol)
    
    # Note: the decelleration test is needed for the Langevin thermostat only. Brownian thermostat is defined
    # over a larger time-step by its concept.
    def check_dissipation_viscous_drag(self, n):
        """
        Check the dissipation relations: the drag terminal velocity tests,
        aka the drift in case of the electrostatics

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """
        system = self.system
        f = np.zeros((2 * n, 3))
        tor = np.zeros((2 * n, 3))
        dip = np.zeros((2 * n, 3))
        tmp_axis = np.zeros((2 * n, 3))
        tol = 1E-10
        if "EXTERNAL_FORCES" in espressomd.features():
            for k in range(2):
                for i in range(n):
                    ind = i + k * n
                    # Just some random forces
                    f[ind,:] = uniform(-0.5, 500., 3)
                    system.part[ind].ext_force = f[ind,:]
                    if "ROTATION" in espressomd.features():
                        # Just some random torques
                        tor[ind,:] = uniform(-0.5, 500., 3)
                        system.part[ind].ext_torque = tor[ind,:]
                        # Let's set the dipole perpendicular to the torque
                        if "DIPOLES" in espressomd.features():
                            # 2 types of particles correspond to 2 different
                            # perpendicular vectors
                            if ind % 2 == 0:
                                dip[ind,:] = 0.0, tor[ind, 2], -tor[ind, 1]
                            else:
                                dip[ind,:] = -tor[ind, 2], 0.0, tor[ind, 0]
                            system.part[ind].dip = dip[ind,:]
                            # 3rd dynamic axis
                            tmp_axis[ind,:] = np.cross(tor[ind,:], dip[ind,:]) \
                                / (np.linalg.norm(tor[ind]) * np.linalg.norm(dip[ind]))
            # Small number of steps is enough for the terminal velocity within the BD by its definition.
            # A simulation of the original saturation of the velocity.
            system.integrator.run(7)
            system.time = 0.0
            for i in range(n):
                for k in range(2):
                    ind = i + k * n
                    system.part[ind].pos = np.zeros((3))
                    if "DIPOLES" in espressomd.features():
                        system.part[ind].dip = dip[ind,:]
            for step in range(3):
                # Small number of steps
                system.integrator.run(2)
                for k in range(2):
                    ind = i + k * n
                    for j in range(3):
                        # Eq. (14.34) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010)
                        # First (deterministic) term of the eq. (14.34) of Schlick2010 taking into account eq. (14.35).
                        self.assertLess(
                            abs(system.part[ind].v[j] - f[ind, j] / \
                                self.gamma_tran_p_validate[k, j]), tol)
                        # Second (deterministic) term of the Eq. (14.39) of Schlick2010.
                        self.assertLess(
                            abs(system.part[ind].pos[j] - \
                                system.time * f[ind, j] / self.gamma_tran_p_validate[k, j]), tol)
                        # Same, a rotational analogy.
                        if "ROTATION" in espressomd.features():
                            self.assertLess(abs(
                                system.part[ind].omega_lab[j] - tor[ind, j] \
                                    / self.gamma_rot_p_validate[k, j]), tol)
                    if "ROTATION" in espressomd.features() and "DIPOLES" in espressomd.features():
                        # Same, a rotational analogy. One is implemented using a simple linear algebra;
                        # the polar angles with a sign control just for a correct inverse trigonometric functions application.
                        cos_alpha = np.dot(dip[ind,:], system.part[ind].dip[:]) / \
                            (np.linalg.norm(dip[ind,:]) * system.part[ind].dipm)
                        # Isoptropic particle for the BD. Single gamma equals to other components
                        cos_alpha_test = np.cos(system.time * np.linalg.norm(tor[ind,:]) / \
                            self.gamma_rot_p_validate[k, 0])
                        # The sign instead of sin calc additionally (equivalent approach)
                        sgn = np.sign(np.dot(system.part[ind].dip[:], tmp_axis[ind,:]))
                        sgn_test = np.sign(np.sin(system.time * np.linalg.norm(tor[ind,:]) / \
                            self.gamma_rot_p_validate[k, 0]))

                        self.assertLess(abs(cos_alpha - cos_alpha_test), tol)
                        self.assertEqual(sgn, sgn_test)

    def check_fluctuation_dissipation(self, n, therm_steps, loops):
        """
        Check the fluctuation-dissipation relations: thermalization
        and diffusion properties.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.
        therm_steps : :obj:`int`
            Number of thermalization steps.
        loops : :obj:`int`
            Number of iterations in the sampling.

        """
        system = self.system
        # The thermalization and diffusion test
        # Checks if every degree of freedom has 1/2 kT of energy, even when
        # mass and inertia tensor are active
        # Check the factual translational diffusion.
        #
        # matrices: [2 types of particless] x [3 dimensions X Y Z]
        # velocity^2, omega^2, position^2
        v2 = np.zeros((2, 3))
        o2 = np.zeros((2, 3))
        dr2 = np.zeros((2, 3))
        # Variance to compare with:
        sigma2_tr = np.zeros((2))
        # Comparable variance:
        dr_norm = np.zeros((2))

        pos0 = np.zeros((2 * n, 3))
        for p in range(n):
            for k in range(2):
                ind = p + k * n
                pos0[ind,:] = system.part[ind].pos
        dt0 = self.mass / self.gamma_tran_p_validate

        system.integrator.run(therm_steps)

        int_steps = 50
        for i in range(loops):
            system.integrator.run(int_steps)
            # Get kinetic energy in each degree of freedom for all particles
            for p in range(n):
                for k in range(2):
                    ind = p + k * n
                    v = system.part[ind].v
                    if "ROTATION" in espressomd.features():
                        o = system.part[ind].omega_body
                        o2[k,:] = o2[k,:] + np.power(o[:], 2)
                    pos = system.part[ind].pos
                    v2[k,:] = v2[k,:] + np.power(v[:], 2)
                    dr2[k,:] = np.power((pos[:] - pos0[ind,:]), 2)
                    dt = (int_steps * (i + 1) + therm_steps) * \
                        system.time_step
                    # translational diffusion variance: after a closed-form
                    # integration of the Langevin EOM;
                    # ref. the eq. (10.2.26) N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010)
                    # after simple transformations and the dimensional model
                    # matching (cf. eq. (10.1.1) there):
                    sigma2_tr[k] = 0.0
                    for j in range(3):
                        sigma2_tr[k] += \
                            self.D_tran_p_validate[k, j] * \
                            (2.0 * dt + dt0[k, j] * \
                            (- 3.0 + 4.0 * \
                            math.exp(- dt / dt0[k, j]) - \
                            math.exp(- 2.0 * \
                            dt / dt0[k, j])))
                    dr_norm[k] += (sum(dr2[k,:]) -
                                   sigma2_tr[k]) / sigma2_tr[k]

        tolerance = 0.15
        Ev = 0.5 * self.mass * v2 / (n * loops)
        Eo = 0.5 * self.J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k,:]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do[k] = sum(Eo[k,:]) / (3.0 * self.halfkT_p_validate[k]) - 1.0
            do_vec[k,:] = Eo[k,:] / self.halfkT_p_validate[k] - 1.0
        dr_norm /= (n * loops)

        for k in range(2):
            self.assertLessEqual(
                abs(
                    dv[k]),
                tolerance,
                msg='Relative deviation in translational energy too large: {0}'.format(
                    dv[k]))
            if "ROTATION" in espressomd.features():
                self.assertLessEqual(
                    abs(
                        do[k]),
                    tolerance,
                    msg='Relative deviation in rotational energy too large: {0}'.format(
                        do[k]))
                self.assertLessEqual(abs(
                    do_vec[k, 0]), tolerance, msg='Relative deviation in rotational energy per the body axis X is too large: {0}'.format(do_vec[k, 0]))
                self.assertLessEqual(abs(
                    do_vec[k, 1]), tolerance, msg='Relative deviation in rotational energy per the body axis Y is too large: {0}'.format(do_vec[k, 1]))
                self.assertLessEqual(abs(
                    do_vec[k, 2]), tolerance, msg='Relative deviation in rotational energy per the body axis Z is too large: {0}'.format(do_vec[k, 2]))
            self.assertLessEqual(
                abs(
                    dr_norm[k]),
                tolerance,
                msg='Relative deviation in translational diffusion is too large: {0}'.format(
                    dr_norm[k]))

    def set_particle_specific_gamma(self, n):
        """
        Set the particle-specific gamma.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            # for the expected metrics calc
            self.gamma_tran_p_validate[k,:] = self.gamma_tran_p[k,:]
            self.gamma_rot_p_validate[k,:] = self.gamma_rot_p[k,:]
            # init
            for i in range(n):
                ind = i + k * n
                self.system.part[ind].gamma = self.gamma_tran_p[k,:]
                if "ROTATION" in espressomd.features():
                    self.system.part[ind].gamma_rot = self.gamma_rot_p[k,:]

    def set_particle_specific_temperature(self, n):
        """
        Set the particle-specific temperature.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles of the each type. There are 2 types.

        """

        for k in range(2):
            # expected
            self.halfkT_p_validate[k] = self.kT_p[k] / 2.0
            # init
            for i in range(n):
                ind = i + k * n
                self.system.part[ind].temp = self.kT_p[k]

    def set_diffusivity_tran(self):
        """
        Set the translational diffusivity to validate further.

        """

        for k in range(2):
            # Translational diffusivity for a validation
            self.D_tran_p_validate[k,:] = 2.0 * \
                self.halfkT_p_validate[k] / self.gamma_tran_p_validate[k,:]

    # Test case 0.0.0:
    # no particle specific values / dissipation only / LD only.
    # No meaning for the simple BD propagation cause
    # it has no inertial features.
    def test_case_000(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        # Actual integration and validation run
        self.check_dissipation(n)

    # Test case 0.0.1:
    # no particle specific values / dissipation viscous drag only / BD only.
    # LD will require too much computational time
    # (one is tested offline though).
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_001(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 0.0.2:
    # no particle specific values / dissipation viscous drag only / EB only.
    # LD will require too much computational time
    # (one is tested offline though).
    if "ERMAK_BUCKHOLZ" in espressomd.features():
        def test_case_002(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 0.1: no particle specific values / fluctuation & dissipation
    # Same particle and thermostat parameters for LD and BD/EB are required in order
    # to test the BD/EB consistency by means of NVT-ensemble.
    # Less number of steps and other specific integration parameters of BD/EB
    # reflect its temporal scale advances.
    def test_case_01(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 500
        therm_steps = 1000
        loops = 10
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            # Large time-step is OK for BD.
            system.time_step = 10.0
            # Less number of loops are needed in case of BD because the velocity
            # distribution is already as required. It is not a result of a real dynamics.
            loops = 8
            # The BD does not require so the warmup. Only 1 step is enough.
            # More steps are taken just to be sure that they will not lead
            # to wrong results.
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "ERMAK_BUCKHOLZ" in espressomd.features():
            self.set_initial_cond()
            # Large time-step is OK for EB.
            system.time_step = 10.0
            # Less number of loops are needed in case of BD because the velocity
            # distribution is already as required. It is not a result of a real dynamics.
            loops = 8
            # The EB does not require so the warmup. Only 1 step is enough.
            # More steps are taken just to be sure that they will not lead
            # to wrong results.
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 1.0.0: particle specific gamma but not temperature / dissipation
    # only / LD only
    def test_case_100(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_gamma(n)
        # Actual integration and validation run
        self.check_dissipation(n)

    # Test case 1.0.1: particle specific gamma but not temperature /
    # dissipation viscous drag only / BD only.
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_101(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_gamma(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 1.0.2: particle specific gamma but not temperature /
    # dissipation viscous drag only / EB only.
    if "ERMAK_BUCKHOLZ" in espressomd.features():
        def test_case_102(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_gamma(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 1.1: particle specific gamma but not temperature / fluctuation
    # & dissipation / LD and BD/EB
    def test_case_11(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 500
        therm_steps = 1000
        loops = 10
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_gamma(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "ERMAK_BUCKHOLZ" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 2.0.0: particle specific temperature but not gamma / dissipation
    # only / LD only
    def test_case_200(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation(n)

    # Test case 2.0.1: particle specific temperature but not gamma / dissipation
    # viscous drag only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_201(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 2.0.2: particle specific temperature but not gamma / dissipation
    # viscous drag only / EB only
    if "ERMAK_BUCKHOLZ" in espressomd.features():
        def test_case_202(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 2.1: particle specific temperature but not gamma / fluctuation
    # & dissipation / LD and BD/EB
    def test_case_21(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 500
        therm_steps = 1000
        loops = 10
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "ERMAK_BUCKHOLZ" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.turn_off()
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 3.0.0: both particle specific gamma and temperature /
    # dissipation only / LD only
    def test_case_300(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        # Actual integration and validation run
        self.check_dissipation(n)

    # Test case 3.0.1: both particle specific gamma and temperature /
    # dissipation viscous drag only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_301(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_gamma(n)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 3.0.2: both particle specific gamma and temperature /
    # dissipation viscous drag only / EB only
    if "ERMAK_BUCKHOLZ" in espressomd.features():
        def test_case_302(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.set_particle_specific_gamma(n)
            self.set_particle_specific_temperature(n)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 3.1: both particle specific gamma and temperature /
    # fluctuation & dissipation / LD and BD/EB
    def test_case_31(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 500
        therm_steps = 1000
        loops = 10
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        self.set_particle_specific_gamma(n)
        self.set_particle_specific_temperature(n)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "ERMAK_BUCKHOLZ" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_eb(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

    # Test case 4.0.0: no particle specific values / rotational specific global
    # thermostat / dissipation only / LD only
    def test_case_400(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 1
        self.dissipation_param_setup(n)
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        self.system.thermostat.turn_off()
        system.thermostat.set_langevin(
            kT=self.kT,
            gamma=self.gamma_global,
            gamma_rotation=self.gamma_global_rot, seed=42)
        # Actual integration and validation run
        self.check_dissipation(n)

    # Test case 4.0.1: no particle specific values / rotational specific global
    # thermostat / dissipation only / BD only
    if "BROWNIAN_DYNAMICS" in espressomd.features():
        def test_case_401(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults_rot_differ()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_brownian(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot,
                seed=42)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 4.0.2: no particle specific values / rotational specific global
    # thermostat / dissipation only / EB only
    if "ERMAK_BUCKHOLZ" in espressomd.features():
        def test_case_402(self):
            system = self.system
            # Each of 2 kind of particles will be represented by n instances:
            n = 1
            self.dissipation_viscous_drag_setup_bd()
            self.set_langevin_global_defaults_rot_differ()
            # The test case-specific thermostat and per-particle parameters
            system.thermostat.set_eb(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot,
                seed=42)
            # Actual integration and validation run
            self.check_dissipation_viscous_drag(n)

    # Test case 4.1: no particle specific values / rotational specific global
    # thermostat / fluctuation & dissipation / LD and BD/EB
    def test_case_41(self):
        system = self.system
        # Each of 2 kind of particles will be represented by n instances:
        n = 500
        therm_steps = 1000
        loops = 10
        self.fluctuation_dissipation_param_setup(n)
        self.set_langevin_global_defaults_rot_differ()
        # The test case-specific thermostat and per-particle parameters
        system.thermostat.set_langevin(
            kT=self.kT,
            gamma=self.gamma_global,
            gamma_rotation=self.gamma_global_rot, seed=42)
        self.set_diffusivity_tran()
        # Actual integration and validation run
        self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_brownian(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot,
                seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)
        if "ERMAK_BUCKHOLZ" in espressomd.features():
            self.set_initial_cond()
            system.time_step = 10.0
            loops = 8
            therm_steps = 2
            # The test case-specific thermostat
            system.thermostat.turn_off()
            system.thermostat.set_eb(
                kT=self.kT,
                gamma=self.gamma_global,
                gamma_rotation=self.gamma_global_rot,
                seed=42)
            # Actual integration and validation run
            self.check_fluctuation_dissipation(n, therm_steps, loops)

if __name__ == '__main__':
    ut.main()
