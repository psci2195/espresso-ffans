#ifndef INTEGRATORS_VELOCITY_VERLET_HPP
#define INTEGRATORS_VELOCITY_VERLET_HPP

#include "ParticleRange.hpp"
#include "config.hpp"
#include "particle_data.hpp"
#include "rotation.hpp"

/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
 *  v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+0.5 \Delta t) \f]
 */
inline void velocity_verlet_propagate_vel_pos(ParticleRange &particles) {

  auto const skin2 = Utils::sqr(0.5 * skin);
  for (auto &p : particles) {
#ifdef ROTATION
    propagate_omega_quat_particle(p);
#endif

// Don't propagate translational degrees of freedom of vs
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif

#ifdef GRONBECH_JENSEN_FARAGO
    // the lab system for sure:
    Utils::Vector3d v_0_buf;
    Utils::Vector3d p_0_buf;
    if (integ_switch == INTEG_METHOD_GRONBECH_J_FARAGO) {
      // save initial conditions
      v_0_buf = p.m.v;
      p_0_buf = p.r.p;
    }
#endif // GRONBECH_JENSEN_FARAGO

    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
      {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
         * v(t+0.5*dt) */
        p.r.p[j] += time_step * p.m.v[j];
      }
    }

#ifdef GRONBECH_JENSEN_FARAGO
    //Utils::Vector3d f_0_buf;
    if (integ_switch == INTEG_METHOD_GRONBECH_J_FARAGO) {
      Thermostat::GammaType langevin_pref_friction_buf;
      Thermostat::GammaType langevin_pref_noise_buf;
      bool aniso_flag;
      aniso_flag = false;
      friction_thermo_langevin_pref(p, langevin_pref_friction_buf,
                                    langevin_pref_noise_buf, aniso_flag);
      // possible body system of references:
      Utils::Vector3d velocity_buf;
      Utils::Vector3d position_buf;
      Utils::Vector3d velocity_0_buf;
      Utils::Vector3d position_0_buf;

      p.m.noise_saved_0 = p.m.noise_saved;
      //f_0_buf = p.f.f;
#ifdef PARTICLE_ANISOTROPY
      if (aniso_flag) {
        velocity_buf = convert_vector_space_to_body(p, p.m.v);
        position_buf = convert_vector_space_to_body(p, p.r.p);
        velocity_0_buf = convert_vector_space_to_body(p, v_0_buf);
        position_0_buf = convert_vector_space_to_body(p, p_0_buf);
        //f_0_buf = convert_vector_space_to_body(p, p.f.f);
      } else
#endif // PARTICLE_ANISOTROPY
      {
        // the vectors are set within the space system by default
        velocity_buf = p.m.v;
        position_buf = p.r.p;
        velocity_0_buf = v_0_buf;
        position_0_buf = p_0_buf;
      }

      for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
        if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
        {
          double b = 1. /(1. - langevin_pref_friction_buf[j] * 
                          time_step / (2. * p.p.mass));
          // get read of the friction part of the velocity:
          velocity_buf[j] -= 0.5 * time_step * langevin_pref_friction_buf[j] *
                              velocity_0_buf[j] / p.p.mass; 
          //f_0_buf[j] -= langevin_pref_friction_buf[j] * velocity_0_buf[j];
          // cleanup the double positional update from the Step 1 VV:
          position_buf[j] -= velocity_0_buf[j] * time_step;
          // now, the main propagation. We assume the VV Step 1 changes above made!
          position_buf[j] += b * velocity_buf[j] * time_step;
          velocity_buf[j] += (langevin_pref_friction_buf[j] / p.p.mass) *
                              (position_buf[j] - position_0_buf[j]);
        }
      }

    // transform the phase space coord. back to the lab system if needed
    if (aniso_flag) {
        p.m.v = convert_vector_body_to_space(p, velocity_buf);
        p.r.p = convert_vector_body_to_space(p, position_buf);
      } else {
        p.m.v = velocity_buf;
        p.r.p = position_buf;
      }
    }
#endif // GRONBECH_JENSEN_FARAGO

    /* Verlet criterion check*/
    if (Utils::sqr(p.r.p[0] - p.l.p_old[0]) +
            Utils::sqr(p.r.p[1] - p.l.p_old[1]) +
            Utils::sqr(p.r.p[2] - p.l.p_old[2]) >
        skin2)
      set_resort_particles(Cells::RESORT_LOCAL);
  }
}

/** Final integration step of the Velocity Verlet integrator
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f]
 */
inline void
velocity_verlet_propagate_vel_final(const ParticleRange &particles) {

  for (auto &p : particles) {
#ifdef VIRTUAL_SITES
    // Virtual sites are not propagated during integration
    if (p.p.is_virtual)
      continue;
#endif

#ifdef GRONBECH_JENSEN_FARAGO
    // the lab system for sure:
    Utils::Vector3d v_0_buf;
    Utils::Vector3d noise_0_buf;
    if (integ_switch == INTEG_METHOD_GRONBECH_J_FARAGO) {
      // save initial conditions (for this step!)
      v_0_buf = p.m.v;
      noise_0_buf = p.m.noise_saved;
    }
#endif // GRONBECH_JENSEN_FARAGO

    for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & COORD_FIXED(j))) {
#endif
        /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
        p.m.v[j] += 0.5 * time_step * p.f.f[j] / p.p.mass;
#ifdef EXTERNAL_FORCES
      }
#endif
    }

#ifdef GRONBECH_JENSEN_FARAGO
    //Utils::Vector3d f_0_buf;
    if (integ_switch == INTEG_METHOD_GRONBECH_J_FARAGO) {
      Thermostat::GammaType langevin_pref_friction_buf;
      Thermostat::GammaType langevin_pref_noise_buf;
      bool aniso_flag;
      aniso_flag = false;
      friction_thermo_langevin_pref(p, langevin_pref_friction_buf,
                                    langevin_pref_noise_buf, aniso_flag);
      // possible body system of references:
      Utils::Vector3d velocity_buf;
      Utils::Vector3d velocity_0_buf;
      Utils::Vector3d noise_0_body_buf;
      Utils::Vector3d noise_body_buf;
      //f_0_buf = p.f.f;
#ifdef PARTICLE_ANISOTROPY
      if (aniso_flag) {
        velocity_buf = convert_vector_space_to_body(p, p.m.v);
        velocity_0_buf = convert_vector_space_to_body(p, v_0_buf);
        noise_body_buf = convert_vector_space_to_body(p, p.m.noise_saved);
        // the noise stored from the Step 1 of the VV:
        noise_0_body_buf = convert_vector_space_to_body(p, p.m.noise_saved_0);
        //f_0_buf = convert_vector_space_to_body(p, p.f.f);
      } else
#endif // PARTICLE_ANISOTROPY
      {
        // the vectors are set within the space system by default
        velocity_buf = p.m.v;
        velocity_0_buf = v_0_buf;
        noise_body_buf = p.m.noise_saved;
        noise_0_body_buf = p.m.noise_saved_0;
      }

      for (int j = 0; j < 3; j++) {
#ifdef EXTERNAL_FORCES
        if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
        {
          // get read of the friction part of the velocity:
          velocity_buf[j] -= 0.5 * time_step * langevin_pref_friction_buf[j] *
                              velocity_0_buf[j] / p.p.mass; 
          // same for the noise:
          velocity_buf[j] -= noise_body_buf[j] * langevin_pref_noise_buf[j] *
                             time_step_half / p.p.mass;
          //f_0_buf[j] -= langevin_pref_friction_buf[j] * velocity_0_buf[j];
          // now, the main propagation. We assume the VV Step 2 changes above made!
          velocity_buf[j] += noise_0_body_buf[j] * langevin_pref_noise_buf[j] *
                             time_step_half / p.p.mass;
        }
      }

    // transform the phase space coord. back to the lab system if needed
    if (aniso_flag) {
        p.m.v = convert_vector_body_to_space(p, velocity_buf);
      } else {
        p.m.v = velocity_buf;
      }
    }
#endif // GRONBECH_JENSEN_FARAGO

  }
}

inline void velocity_verlet_step_1(ParticleRange &particles) {
  velocity_verlet_propagate_vel_pos(particles);
  sim_time += time_step;
}

inline void velocity_verlet_step_2(const ParticleRange &particles) {
  velocity_verlet_propagate_vel_final(particles);
#ifdef ROTATION
  convert_torques_propagate_omega(particles);
#endif
}

#endif
