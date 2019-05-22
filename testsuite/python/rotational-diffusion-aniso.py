from __future__ import print_function
import math
import numpy as np
from numpy.random import seed, uniform
import unittest as ut
import random
import espressomd
import tests_common


@ut.skipIf(not espressomd.has_features(["ROTATION",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "DIPOLES"]),
           "Features not available, skipping test!")
class RotDiffAniso(ut.TestCase):
    longMessage = True
    round_error_prec = 1E-14
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))

    # Particle properties
    J = np.zeros((3))

    @classmethod
    def setUpClass(cls):
        # Handle a random generator seeding
        #rnd_gen = random.SystemRandom()
        #seed1 = int(200 * rnd_gen.random())
        #seed1 = 42
        #np.random.seed(seed1)
        #seed2 = int(200 * rnd_gen.random())
        #seed2 = 42
        # The Espresso system configuration
        #cls.system.seed = [s * seed2 for s in range(cls.system.cell_system.get_state()["n_nodes"])]
        cls.system.set_random_state_PRNG()
        np.random.seed(seed=cls.system.seed)
        cls.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        cls.system.cell_system.skin = 5.0

    def setUp(self):
        self.system.time = 0.0
        self.system.part.clear()

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

    def add_particles_setup(self, n):
        """
        Adding particles according to the
        previously set parameters.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles.

        """

        for ind in range(n):
            part_pos = np.random.random(3) * self.system.box_l[0]
            self.system.part.add(rotation=(1, 1, 1), id=ind, rinertia=self.J,
                                 pos=part_pos)
            if "ROTATION" in espressomd.features():
                self.system.part[ind].omega_body = 0.0, 0.0, 0.0

    def set_anisotropic_param(self):
        """
        Select parameters for anisotropic particles.

        Parameters
        ----------

        """

        # NVT thermostat
        # Note: here & hereinafter specific variations in the random parameter ranges are related to
        # the test execution duration to achieve the required statistical averages faster.
        # The friction gamma_global should be large enough in order to have the small enough D ~ kT / gamma and
        # to observe the details of the original rotational diffusion: the Perrin1936 (see the reference below) tests are visible
        # only when the diffusive rotation is ~pi due to the exponential temporal dependencies (see the equations referred in the check_rot_diffusion()).
        # Also, t0 ~ J / gamma should be small enough
        # in order to remove the small-time-scale diffusion effects which do not fit the Perrin1936's
        # tests which are based on the partial differential equation (eq. (68), Perrin1934) leading only to the simple
        # classical Einstein-Smoluchowski equations of the diffusion
        # in a contrast of the eq. (10.2.26) [N. Pottier,
        # https://doi.org/10.1007/s10955-010-0114-6 (2010)].
        self.gamma_global = 1E2 * uniform(0.35, 1.05, (3))

        # Particles' properties
        # As far as the problem characteristic time is t0 ~ J / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the moment of inertia higher than
        # some minimal value.
        # Also, it is expected to test the large enough J.
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time: the in silico time should clock over the
        # t0.
        self.J = uniform(1.5, 16.5, (3))

    def set_isotropic_param(self):
        """
        Select parameters for isotropic particles.

        Parameters
        ----------

        """

        # NVT thermostat
        # see the comments in set_anisotropic_param()
        self.gamma_global[0] = 1E2 * uniform(0.35, 1.05, (1))
        self.gamma_global[1] = self.gamma_global[0]
        self.gamma_global[2] = self.gamma_global[0]
        # Particles' properties
        # see the comments in set_anisotropic_param()
        self.J[0] = uniform(1.5, 16.5, (1))
        self.J[1] = self.J[0]
        self.J[2] = self.J[0]

    def rot_diffusion_param_setup(self):
        """
        Setup the parameters for the rotational diffusion
        test check_rot_diffusion().

        Parameters
        ----------

        """

        # Time
        # The time step should be less than t0 ~ mass / gamma
        # self.system.time_step = 3E-3
        self.system.time_step = 3E-4

        # Space
        box = 10.0
        self.system.box_l = box, box, box
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.system.periodicity = 0, 0, 0

        # NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = uniform(1.5, 6.5)

    def check_rot_diffusion(self, n):
        """
        The rotational diffusion tests based on the reference work
        [Perrin, F. (1936) Journal de Physique et Le Radium, 7(1), 1-11. https://doi.org/10.1051/jphysrad:01936007010100]
        with a theoretical background of
        [Perrin, F. (1934) Journal de Physique et Le Radium, 5(10), 497-511. https://doi.org/10.1051/jphysrad:01934005010049700]

        Parameters
        ----------
        n : :obj:`int`
            Number of particles.

        """
        # Global diffusivity tensor in the body frame:
        D = self.kT / self.gamma_global
        #dt0 = self.J / self.gamma_global

        # Thermalizing...
        #therm_steps = 100
        therm_steps = int(1E3)
        self.system.integrator.run(therm_steps)

        # Measuring...
        # Set the initial conditions according to the [Perrin1936], p.3.
        # The body angular velocity is rotated now, but there is
        # only the thermal velocity, hence, this does not impact the test and
        # its physical context.
        for ind in range(n):
            self.system.part[ind].quat = 1.0, 0.0, 0.0, 0.0
        # Average direction cosines
        # Diagonal ones:
        dcosjj_validate = np.zeros((3))
        dcosjj_dev = np.zeros((3))
        # same to the power of 2
        dcosjj2_validate = np.zeros((3))
        dcosjj2_dev = np.zeros((3))
        # The non-diagonal elements for 2 different tests: negative ("nn") and
        # positive ("pp") ones.
        dcosijpp_validate = np.ones((3, 3))
        dcosijpp_dev = np.zeros((3, 3))
        dcosijnn_validate = np.ones((3, 3))
        dcosijnn_dev = np.zeros((3, 3))
        # The non-diagonal elements to the power of 2
        dcosij2_validate = np.ones((3, 3))
        dcosij2_dev = np.zeros((3, 3))

        self.system.time = 0.0
        int_steps = 20
        loops = 100
        print('\n t | j | i | dcosjj | dcosjj2 | dcosijpp | dcosijnn | dcosij2')
        print('\n ============================================================')
        for step in range(loops):
            self.system.integrator.run(steps=int_steps)
            dcosjj = np.zeros((3))
            dcosjj2 = np.zeros((3))
            dcosijpp = np.zeros((3, 3))
            dcosijnn = np.zeros((3, 3))
            dcosij2 = np.zeros((3, 3))
            for ind in range(n):
                # Just a direction cosines functions averaging..
                dir_cos = tests_common.rotation_matrix_quat(self.system, ind)
                for j in range(3):
                    # the LHS of eq. (23) [Perrin1936].
                    dcosjj[j] += dir_cos[j, j]
                    # the LHS of eq. (32) [Perrin1936].
                    dcosjj2[j] += dir_cos[j, j]**2.0
                    for i in range(3):
                        if i != j:
                            # the LHS of eq. (24) [Perrin1936].
                            dcosijpp[
                                i,
                                j] += dir_cos[
                                    i,
                                    i] * dir_cos[
                                        j,
                                        j] + dir_cos[
                                            i,
                                            j] * dir_cos[
                                                j,
                                                i]
                            # the LHS of eq. (25) [Perrin1936].
                            dcosijnn[
                                i,
                                j] += dir_cos[
                                    i,
                                    i] * dir_cos[
                                        j,
                                        j] - dir_cos[
                                            i,
                                            j] * dir_cos[
                                                j,
                                                i]
                            # the LHS of eq. (33) [Perrin1936].
                            dcosij2[i, j] += dir_cos[i, j]**2.0
            dcosjj /= n
            dcosjj2 /= n
            dcosijpp /= n
            dcosijnn /= n
            dcosij2 /= n

            # Actual comparison.

            tolerance = 0.12
            # Too small values of the direction cosines are out of interest
            # compare to 0..1 range.
            min_value = 0.14

            # Eq. (23) [Perrin1936].
            dcosjj_validate[0] = np.exp(-(D[1] + D[2]) * self.system.time)
            dcosjj_validate[1] = np.exp(-(D[0] + D[2]) * self.system.time)
            dcosjj_validate[2] = np.exp(-(D[0] + D[1]) * self.system.time)
            dcosjj_dev = np.absolute(
                dcosjj - dcosjj_validate) / dcosjj_validate
            for j in range(3):
                if np.absolute(dcosjj_validate[j]) < min_value:
                    dcosjj_dev[j] = 0.0

            # Eq. (24) [Perrin1936].
            dcosijpp_validate[0, 1] = np.exp(
                -(4 * D[2] + D[1] + D[0]) * self.system.time)
            dcosijpp_validate[1, 0] = np.exp(
                -(4 * D[2] + D[1] + D[0]) * self.system.time)
            dcosijpp_validate[0, 2] = np.exp(
                -(4 * D[1] + D[2] + D[0]) * self.system.time)
            dcosijpp_validate[2, 0] = np.exp(
                -(4 * D[1] + D[2] + D[0]) * self.system.time)
            dcosijpp_validate[1, 2] = np.exp(
                -(4 * D[0] + D[2] + D[1]) * self.system.time)
            dcosijpp_validate[2, 1] = np.exp(
                -(4 * D[0] + D[2] + D[1]) * self.system.time)
            dcosijpp_dev = np.absolute(
                dcosijpp - dcosijpp_validate) / dcosijpp_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijpp_validate[i, j]) < min_value:
                        dcosijpp_dev[i, j] = 0.0

            # Eq. (25) [Perrin1936].
            dcosijnn_validate[0, 1] = np.exp(-(D[1] + D[0]) * self.system.time)
            dcosijnn_validate[1, 0] = np.exp(-(D[1] + D[0]) * self.system.time)
            dcosijnn_validate[0, 2] = np.exp(-(D[2] + D[0]) * self.system.time)
            dcosijnn_validate[2, 0] = np.exp(-(D[2] + D[0]) * self.system.time)
            dcosijnn_validate[1, 2] = np.exp(-(D[2] + D[1]) * self.system.time)
            dcosijnn_validate[2, 1] = np.exp(-(D[2] + D[1]) * self.system.time)
            dcosijnn_dev = np.absolute(
                dcosijnn - dcosijnn_validate) / dcosijnn_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijnn_validate[i, j]) < min_value:
                        dcosijnn_dev[i, j] = 0.0

            # Eq. (30) [Perrin1936].
            D0 = sum(D[:]) / 3.0
            D1D1 = 0.0
            for j in range(3):
                for i in range(3):
                    if i != j:
                        D1D1 += D[i] * D[j]
            D1D1 /= 6.0
            # Technical workaround of a digital arithmetic issue for isotropic particle
            if np.absolute((D0**2 - D1D1)/(D0**2 + D1D1)) < self.round_error_prec:
                D1D1 *= (1.0 - 2.0 * self.round_error_prec)
            # Eq. (32) [Perrin1936].
            dcosjj2_validate = 1. / 3. + (1. / 3.) * (1. + (D - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                + (1. / 3.) * (1. - (D - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosjj2_dev = np.absolute(
                dcosjj2 - dcosjj2_validate) / dcosjj2_validate
            for j in range(3):
                if np.absolute(dcosjj2_validate[j]) < min_value:
                    dcosjj2_dev[j] = 0.0

            # Eq. (33) [Perrin1936].
            dcosij2_validate[0, 1] = 1. / 3. - (1. / 6.) * (1. - (D[2] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[2] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[1, 0] = dcosij2_validate[0, 1]

            dcosij2_validate[0, 2] = 1. / 3. - (1. / 6.) * (1. - (D[1] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[1] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[2, 0] = dcosij2_validate[0, 2]

            dcosij2_validate[1, 2] = 1. / 3. - (1. / 6.) * (1. - (D[0] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[0] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[2, 1] = dcosij2_validate[1, 2]
            dcosij2_dev = np.absolute(
                dcosij2 - dcosij2_validate) / dcosij2_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosij2_validate[i, j]) < min_value:
                        dcosij2_dev[i, j] = 0.0

            for j in range(3):
                self.assertLessEqual(
                    abs(
                        dcosjj_dev[j]),
                    tolerance,
                    msg='Relative deviation dcosjj_dev[{0}] in a rotational diffusion is too large: {1}'.format(
                        j,
                        dcosjj_dev[j]))
                self.assertLessEqual(
                    abs(
                        dcosjj2_dev[j]),
                    tolerance,
                    msg='Relative deviation dcosjj2_dev[{0}] in a rotational diffusion is too large: {1}'.format(
                        j,
                        dcosjj2_dev[j]))
                for i in range(3):
                    print('\n {0},{1},{2},{3},{4},{5},{6},{7}'.format(
                        self.system.time,j,i,
                        dcosjj[j]/dcosjj_validate[j],
                        dcosjj2[j]/dcosjj2_validate[j],
                        dcosijpp[i, j]/dcosijpp_validate[i, j],
                        dcosijnn[i, j]/dcosijnn_validate[i, j],
                        dcosij2[i, j]/dcosij2_validate[i, j]
                        ))
                    if i != j:
                        self.assertLessEqual(
                            abs(
                                dcosijpp_dev[i, j]),
                            tolerance,
                            msg='Relative deviation dcosijpp_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i, j, dcosijpp_dev[i, j]))
                        self.assertLessEqual(
                            abs(
                                dcosijnn_dev[i, j]),
                            tolerance,
                            msg='Relative deviation dcosijnn_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i, j, dcosijnn_dev[i, j]))
                        self.assertLessEqual(
                            abs(
                                dcosij2_dev[i, j]),
                            tolerance,
                            msg='Relative deviation dcosij2_dev[{0},{1}] in a rotational diffusion is too large: {2}'.format(
                                i, j, dcosij2_dev[i, j]))

    def state_print(self, check):
        system = self.system
        print('\n \n', check)
        #print('\n', system.thermostat.get_state())
        thermo_state = system.thermostat.get_state()
        print('\n kT={0}, gamma={1}, type={2}'.format(thermo_state[0]["kT"], thermo_state[0]["gamma"], thermo_state[0]["type"]))
        part = system.part[0]
        print('\n mass={0} rintertia={1}'.format(part.mass, part.rinertia))

    def test_case_00(self):
        n = int(1.E3)
        self.rot_diffusion_param_setup()
        self.set_anisotropic_param()
        self.add_particles_setup(n)
        self.system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.state_print(check = 'LD: check_rot_diffusion (aniso)')
        self.check_rot_diffusion(n)

    def test_case_01(self):
        n = int(1.E3)
        self.rot_diffusion_param_setup()
        self.set_isotropic_param()
        self.add_particles_setup(n)
        self.system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.state_print(check = 'LD: check_rot_diffusion (iso)')
        self.check_rot_diffusion(n)
        self.set_initial_cond()
        self.system.thermostat.turn_off()
        self.system.thermostat.set_brownian(
            kT=self.kT, gamma=self.gamma_global)
        # Actual integration and validation run
        self.state_print(check = 'BD: check_rot_diffusion (iso)')
        self.check_rot_diffusion(n)

if __name__ == '__main__':
    ut.main()
