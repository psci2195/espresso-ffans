#
# Copyright (C) 2013-2019 The ESPResSo project
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
#
from __future__ import print_function
import unittest as ut
import numpy as np
from time import time
import random
import espressomd
from espressomd.interactions import HarmonicBond
from espressomd.accumulators import Correlator
from espressomd.observables import ParticleVelocities, ParticleBodyAngularVelocities, ParticlePositions
from tests_common import single_component_maxwell


@ut.skipIf(espressomd.has_features("THERMOSTAT_IGNORE_NON_VIRTUAL"),
           "Skipped because of THERMOSTAT_IGNORE_NON_VIRTUAL")
class HarmonicOscillatorThermalization(ut.TestCase):

    """Tests that Langevin thermostat applies to the harmonic oscillator
    accordig to the statistical physics."""

    system = None
    hb = None

    @classmethod
    def setUpClass(cls):
        cls.system = espressomd.System(box_l=[1.0E2, 1.0E2, 1.0E2])
        if espressomd.has_features("PARTIAL_PERIODIC"):
            cls.system.periodicity = 0, 0, 0
        # Handle a random generator seeding
        #rnd_gen = random.SystemRandom()
        #seed1 = int(200 * rnd_gen.random())
        #seed1 = 15
        #np.random.seed(seed1)
        #seed2 = int(200 * rnd_gen.random())
        #seed2 = 42
        # The Espresso system configuration
        #cls.system.seed = [s * seed2 for s in range(cls.system.cell_system.get_state()["n_nodes"])]
        cls.system.set_random_state_PRNG()
        np.random.seed(seed=cls.system.seed)
        cls.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        cls.system.cell_system.skin = 5.0

        # Harmonic interaction
        hb_k = 1.
        hb_r_0 = 0.0
        hb_r_cut = 20
        cls.hb = espressomd.interactions.HarmonicBond(
            k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)
        cls.system.bonded_inter.add(cls.hb)

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """check the recorded particle distributions in vel againsta histogram with n_bins bins. Drop velocities outside minmax. Check individual histogram bins up to an accuracy of error_tol agaisnt the analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, normed=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            for j in range(n_bins):
                found = data[j]
                expected = single_component_maxwell(
                    bins[j], bins[j + 1], kT)
                self.assertLessEqual(abs(found - expected), error_tol)

    def test_aa_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertLessEqual(
            abs(single_component_maxwell(-10, 10, 4.) - 1.), 1E-4)

    def global_langevin_run_check(self, N, kT, loops, error_tol):
        """Sampling for the global Langevin parameters test.

        Parameters
        ----------
        N :       :obj:`int`
                  Number of particles.
        kT :      :obj:`float`
                  Temperature.
        loops :   :obj:`int`
                  Number of integration steps to sample.

        """
        system = self.system
        v_stored = np.zeros((N * loops, 3))
        omega_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            system.integrator.run(1)
            v_stored[i * N:(i + 1) * N, :] = system.part[:].v
            if espressomd.has_features("ROTATION"):
                omega_stored[i * N:(i + 1) * N, :] = system.part[:].omega_body

        v_minmax = 5
        bins = 4
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)
        if espressomd.has_features("ROTATION"):
            self.check_velocity_distribution(
                omega_stored, v_minmax, bins, error_tol, kT)

    def harmonic_set_pairs(self, N):
        """Set the harmonic interaction between
        odd and even particles."""
        # Harmonic interactions
        for i in range(N // 2):
            self.system.part[i].add_bond((self.hb, i + N // 2))

    def test_global_langevin(self):
        """Test for global Langevin parameters."""
        # Must be even here
        N = 200
        system = self.system
        system.part.clear()
        system.time_step = 1E-2

        # Place particles
        system.part.add(pos=np.random.random((N, 3)))

        # Enable the harmonic interaction
        self.harmonic_set_pairs(N)

        # Enable rotation if compiled in
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = 1, 1, 1

        kT = 1.
        gamma = 1.
        system.thermostat.set_langevin(kT=kT, gamma=gamma)

        # Warmup
        system.integrator.run(int(4E3))

        self.global_langevin_run_check(N, kT, 400*6, error_tol = 0.050)

        if espressomd.has_features("BROWNIAN_DYNAMICS"):
            # Large time-step is OK for BD.
            system.time_step = 5.E-1
            system.part[:].v = np.zeros((3))
            system.part[:].omega_body = np.zeros((3))
            system.thermostat.turn_off()
            system.thermostat.set_brownian(kT=kT, gamma=gamma)
            # Warmup
            # The BD does not require so the warmup. Only 1 step is enough.
            # More steps are taken just to be sure that they will not lead
            # to wrong results.
            system.integrator.run(3)
            # Less number of loops are needed in case of BD because the
            # velocity distribution is already as required.
            # It is not a result of a real dynamics.
            self.global_langevin_run_check(N, kT, 400*6, error_tol=0.07)
            system.thermostat.turn_off()

    def setup_diff_mass_rinertia(self, p):
        if espressomd.has_features("MASS"):
            # Beard's "k"
            k = 1
            p.mass = k
        if espressomd.has_features("ROTATION"):
            p.rotation = 1, 1, 1
            # Make sure rinertia does not change diff coeff
            if espressomd.has_features("ROTATIONAL_INERTIA"):
                p.rinertia = k, k, k

    def test_diffusion(self):
        """This tests rotational and translational diffusion coeff via
        green-kubo-like integrals for velocities and positions
        according to the harmonic oscillator known features"""
        system = self.system
        system.part.clear()

        kT = 1.
        dt = 0.1
        system.time_step = dt

        # Translational gamma. We cannot test per-component, if rotation is on,
        # because body and space frames become different.
        gamma = 1.

        # Rotational gamma
        gamma_rot_i = 1.
        gamma_rot_a = 1., 1., 1.

        # Particle with global thermostat params
        p_global = system.part.add(pos=[0., 0., 0.])
        system.part.add(pos=[0., 0., 0.], fix=[1, 1, 1])
        # Make sure, mass doesn't change diff coeff
        self.setup_diff_mass_rinertia(p_global)

        # Enable the harmonic interaction
        self.harmonic_set_pairs(2)

        # Thermostat setup
        if espressomd.has_features("ROTATION"):
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                # particle anisotropy and rotation
                system.thermostat.set_langevin(
                    kT=kT, gamma=gamma, gamma_rotation=gamma_rot_a)
            else:
                # Rotation without particle anisotropy
                system.thermostat.set_langevin(
                    kT=kT, gamma=gamma, gamma_rotation=gamma_rot_i)
        else:
            # No rotation
            system.thermostat.set_langevin(kT=kT, gamma=gamma)

        # system.cell_system.skin = 0.4
        system.integrator.run(int(1E4))

        # Correlators
        pos_obs = {}
        vel_obs = {}
        omega_obs = {}
        corr_pos = {}
        corr_vel = {}
        corr_omega = {}

        # linear pos & vel
        pos_obs = ParticlePositions(ids=system.part[:].id)
        vel_obs = ParticleVelocities(ids=system.part[:].id)
        corr_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=20., delta_N=1,
                              corr_operation="componentwise_product", compress1="discard1")
        corr_vel = Correlator(obs1=vel_obs, tau_lin=16, tau_max=20., delta_N=1,
                              corr_operation="componentwise_product", compress1="discard1")
        system.auto_update_accumulators.add(corr_pos)
        system.auto_update_accumulators.add(corr_vel)
        # angular vel
        if espressomd.has_features("ROTATION"):
            omega_obs = ParticleBodyAngularVelocities(ids=system.part[:].id)
            corr_omega = Correlator(
                obs1=omega_obs, tau_lin=16, tau_max=20, delta_N=1,
                                    corr_operation="componentwise_product", compress1="discard1")
            system.auto_update_accumulators.add(corr_omega)

        system.integrator.run(int(6.E6))

        system.auto_update_accumulators.remove(corr_pos)
        corr_pos.finalize()
        system.auto_update_accumulators.remove(corr_vel)
        corr_vel.finalize()
        if espressomd.has_features("ROTATION"):
            system.auto_update_accumulators.remove(corr_omega)
            corr_omega.finalize()

        # Verify diffusion
        # Translation
        # Cast gammas to vector, to make checks independent of
        # PARTICLE_ANISOTROPY
        gamma = np.ones(3) * gamma
        self.verify_diffusion(p_global, corr_vel, kT, gamma)
        self.verify_diffusion_pos(p_global, corr_pos, kT, gamma)

        # Rotation
        if espressomd.has_features("ROTATION"):
            # Decide on effective gamma rotation, since for rotation it is
            # direction dependent
            eff_gamma_rot = None
            if espressomd.has_features("PARTICLE_ANISOTROPY"):
                eff_gamma_rot = gamma_rot_a
            else:
                eff_gamma_rot = gamma_rot_i * np.ones(3)

            #not relevant
            #self.verify_diffusion(p_global, corr_omega, kT, eff_gamma_rot)

    def verify_diffusion(self, p, corr, kT, gamma):
        """Verifify diffusion coeff.

           p: particle, corr: dict containing correltor with particle as key,
           kT=kT, gamma=gamma as 3 component vector.
        """
        c = corr
        # Integral of vacf via Green-Kubo
        # D= int_0^infty <v(t_0)v(t_0+t)> dt     (o 1/3, since we work
        # componentwise)
        i = p.id
        acf = c.result()[:, [0, 2 + 3 * i, 2 + 3 * i + 1, 2 + 3 * i + 2]]
        np.savetxt("acf_vel.dat", acf)

        # Integrate w. trapez rule
        ratio_average = 0.
        for coord in 1, 2, 3:
            I = np.trapz(acf[:, coord], acf[:, 0])
            ratio = I / (kT / gamma[coord - 1])
            self.assertAlmostEqual(ratio, 0., delta=0.07)
            ratio_average += ratio
        ratio_average /= 3.
        print("\n Green-Kubo-velocity: time_step={0} ratio={1}".format(self.system.time_step, ratio_average))

    def verify_diffusion_pos(self, p, corr, kT, gamma):
        c = corr
        i = p.id
        acf = c.result()[:, [0, 2 + 3 * i, 2 + 3 * i + 1, 2 + 3 * i + 2]]
        np.savetxt("acf_pos.dat", acf)

        # Integrate w. trapez rule
        ratio_average = 0.
        for coord in 1, 2, 3:
            I = np.trapz(acf[:, coord], acf[:, 0])
            ratio = I
            self.assertAlmostEqual(ratio, 1., delta=0.07)
            ratio_average += ratio
        ratio_average /= 3.
        print("\n Green-Kubo-like-position: time_step={0} ratio={1}".format(self.system.time_step, ratio_average))

if __name__ == "__main__":
    ut.main()
