#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

import math
import espressomd
from espressomd import visualization
from numpy.random import random
from espressomd.interactions import *
from espressomd.magnetostatics import *
from espressomd.observables import *
from espressomd.analyze import *
from threading import Thread



# Handle for espresso system
es = espressomd.System()

#Constants - math
pi = 3.14159

#[Si] Constants - physics
Na = 6.022E23
mu0 = 4.0 * pi * 1E-7
muB = 9.27400968 * 1E-24
kb = 1.3806488 * 1E-23
Oe_to_Am = 79.577


#[Si]Materials

#Medium
#[Pa*s] Kerosene dynamic viscousity
eta_car = 0.00149
#[J] Hamaker constant
A_h = 4.0E-20
T = 300.0
Hx = 0
Hy = 0
Hz = 2000 #Oersted 

#Nanoparticle
#[m] magnetite unit cell size - a cubic spinel structure with space group Fd3m (above the Verwey temperature) \cite{Fleer1981}
a0 = 0.8397E-9
#[A/m] magnetite saturation magnetization \cite{Kemp2016}
#Ms = 478.0E3
Ms = 412.0E3
#[kg/m^3] magnetite mass density 
rop = 5240.0

# Oleic acid
#[m^-2] Surface density of oleic acid at the 50% coating \cite{Fertman}
N_oa = 1.0E18
#[m] the oleic acid molecule length
delta = 1.97E-9
#[kg] the oleic acid molecule mass
mass_o = 282.47E-3/Na
# Level of the coverage of the particle surface by the oleic acid
k_o = 0.5
# Flocculation parameter \cite{Cebers1987}
k_f = 1.1

#[SI] Ferrofluid
# The volume fraction value of a FF dispersed phase
#phi_v = 1.2E-6

# The large particles fraction \cite{Ivanov1996}
#PHI_l = 7.5E-2
PHI_l = 0.9

#[m] the mean nanoparticle diameter. Hard part: a ferromagnetic core & double nonmagnetic surface
d_small_hard = 11.5E-9
d_small_hydrod = d_small_hard+2.0*delta
d_small_mag_core = d_small_hard-2.0*a0
V_small_hard = pi*pow(d_small_hard,3.0)/6.0
S_small_hard = pi*pow(d_small_hard,2.0)
N_oleic_small = S_small_hard*k_o*N_oa/0.5
V_small_mag_core = pi*pow(d_small_mag_core,3)/6.0
m_small = Ms*V_small_mag_core
M_small = rop*V_small_hard+N_oleic_small*mass_o
Vo_small_eff = (pi/6.0)*(pow(d_small_hydrod,3)-pow(d_small_hard,3))

# Effective mass density of the oleic molecules on the nanoparticle surface
roo_small_eff = N_oleic_small*mass_o/(Vo_small_eff)
I_small = (pi*pow(d_small_hydrod,5)/60.0)*((rop-roo_small_eff)*pow(1-2*delta/d_small_hydrod,5)+roo_small_eff)

#[m] the large nanoparticles diameter \cite{Ivanov1996}
d_l_hard = 21.6E-9
d_l_hydrod = d_l_hard+2*delta
d_l_mag_core = d_l_hard-2*a0
V_l_hard = pi*pow(d_l_hard,3.0)/6.0
S_l_hard = pi*pow(d_l_hard,2.0)
N_oleic_l = S_l_hard*k_o*N_oa/0.5
V_l_mag_core = pi*pow(d_l_mag_core,3)/6.0
m_l = Ms*V_l_mag_core
M_l = rop*V_l_hard+N_oleic_l*mass_o
Vo_l_eff = (pi/6.0)*(pow(d_l_hydrod,3)-pow(d_l_hard,3))

# Effective mass density of the oleic molecules on the nanoparticle surface
roo_l_eff = N_oleic_l*mass_o/(Vo_l_eff)
I_large = (pi*pow(d_l_hydrod,5)/60.0)*((rop-roo_l_eff)*pow(1-2*delta/d_l_hydrod,5)+roo_l_eff)

# Dimensionless model
SIGMA = d_small_hydrod

# Approximate meam mass
M0 = M_small
H0 = pow(kb*T/(4*pi*mu0),0.5)/pow(SIGMA,1.5)
m0 = math.sqrt(4*pi*kb*T*pow(SIGMA,3.0)/mu0)
t0 = math.sqrt(M0*pow(SIGMA,2.0)/(kb*T))
gamma0 = 3*pi*eta_car*pow(SIGMA,2)/math.sqrt(M0*kb*T)

# Computational experiment setup
# time_step should be increased till the stability loss

# => max for pure VV:
# setmd time_step [expr 0.0005]
# => for SEMI_INTEGRATED (AI) method:
# setmd time_step [expr 0.05]

es.time_step = 5E-2
n_part = 2500

#H_ext_Oe = 500.0
H_ext_Oe = 0.0

n_part_small = n_part*(1-PHI_l)
#set box_l [expr pow($n_part*$PI*($PHI_l*pow($d_l_hard/$SIGMA,3.0)+(1-$PHI_l)*pow($d_small_hard/$SIGMA,3.0))/(6.0*$phi_v),1./3.)]
start_lattice_a = 2*d_l_hydrod/SIGMA
buf_l = 10*start_lattice_a
box_l = start_lattice_a*pow(n_part,1/3.0)+2*buf_l

V_carr = pow(box_l,3.0)
box_l_x = box_l*1.0
box_l_y = box_l*1.0
box_l_z = box_l/1.0
es.box_l = [box_l_x,box_l_y,box_l_z]
es.periodicity = [0, 0, 0]
es.cell_system.skin = 0.4
temp = 1
gammat = 1
gammar = 1
es.thermostat.set_langevin(kT = temp, gamma = 1)
es.cell_system.skin = 0
es.cell_system.set_n_square(use_verlet_lists=False)
analyse_step = 1000

#coord_shift = 0.1*start_lattice_a
coord_shift = buf_l
posx = coord_shift
posy = coord_shift
posz = coord_shift
part_dip = np.zeros((3))
l = 15
for i in range(n_part):
    rnd_val_for_p_type = random()
    costheta = 2* random() - 1
    sintheta = math.sin(math.acos(costheta))
    phi = 2.0*pi*random()
    posx = posx+start_lattice_a
    if posx > (box_l_x-buf_l):
        posx = coord_shift
        posy = posy+start_lattice_a
        if posy > (box_l_y-buf_l):
            posy = coord_shift
            posz = posz+start_lattice_a
            if posz > (box_l_z-buf_l):
                    print("Not enough box_l for particles!")

    if rnd_val_for_p_type < PHI_l:
        part_pos = np.array(random(3)) * l
        costheta = 2 * random() - 1
        sintheta = np.sin(np.arcsin(costheta))
        #dx = part_dip[0], dy = part_dip[1] dz = part_dip[2]
        phi = 2 * np.pi * random()
        dipole_modulus = m_l / m0
        part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
        part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
        part_dip[2] = costheta * dipole_modulus
        massval_large = M_l/M0
        I_val_large =I_large/(M0*pow(SIGMA,2))
        gamma_tr_val_large = gamma0*d_l_hydrod/SIGMA
        gamma_rot_val_large = gamma0*pow(d_l_hydrod/SIGMA,3.0)/3.0
        es.part.add(type = 1, id = i, mass = massval_large, rinertia = I_val_large, pos = [posx, posy, posz],dip = part_dip, gamma = gamma_tr_val_large, gamma_rot = gamma_rot_val_large)
    else:
        part_pos = np.array(random(3)) * l
        costheta = 2 * random() - 1
        sintheta = np.sin(np.arcsin(costheta))
        #dx = part_dip[0], dy = part_dip[1] dz = part_dip[2]
        phi = 2 * np.pi * random()
        dipole_modulus = m_small / m0
        part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
        part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
        part_dip[2] = costheta * dipole_modulus
        massval_small = 1.0
        I_val_small = I_small/(M0*pow(SIGMA,2))
        gamma_tr_val_small = gamma0*1.0
        gamma_rot_val_small = gamma0*1.0/3.0
        es.part.add(type = 0, id = i, mass = massval_small, rinertia = I_val_small, pos = [posx, posy, posz],dip = part_dip, gamma = gamma_tr_val_small, gamma_rot = gamma_rot_val_small)

print( "=====================")
print( "Temporal parameters")
print( "=====================")
print( "Small nanoparticles")
print 'd_small_hydrod={0} m'.format(d_small_hydrod)
print 'tau_v_tran= {0} arb. units'.format(massval_small/gamma_tr_val_small)
print 'tau_v_rot={0} arb. units'.format(I_val_small/gamma_rot_val_small)
print( "")
print( "Large nanoparticles")
print  'd_l_hydrod={0} m'.format(d_l_hydrod)
print  'tau_v_tran={0} arb. units'.format(massval_large/gamma_tr_val_large)
print  'tau_v_rot={0} arb. units'.format(I_val_large/gamma_rot_val_large)

print( "[part 10]")

sig = 1.0
cut = 1.12246*sig
eps = 1.5
shift = 0.25
#es.inter.add(0, 0, lennard-jones, eps, sig, cut, shift, 0)
#es.non_bonded_inter[0, 0].lennard_jones.set_params(
 #               epsilon=eps, sigma=sig,
  #              cutoff=cut, shift=shift)

sig = (d_l_hydrod/SIGMA) 
cut = 1.12246*sig
eps = 1.5
shift = 0.25
#es.inter.add(1, 1, lennard-jones, eps, sig, cut, shift, 0)
es.non_bonded_inter[1, 1].lennard_jones.set_params(
                epsilon=eps, sigma=sig,
                cutoff=cut, shift=shift)

sig = ((0.5*d_l_hydrod+0.5*SIGMA)/SIGMA)
cut = 1.12246*sig
eps = 1.5
shift = 0.25
#es.inter.add(0, 1, lennard-jones, eps, sig, cut, shift, 0)
es.non_bonded_inter[0, 1].lennard_jones.set_params(
                epsilon=eps, sigma=sig,
                cutoff=cut, shift=shift)
if "ROTATION" in espressomd.features():
    deg_free = 6
else:
    deg_free = 3 
#prepare_vmd_connection("vmdfile", 10000) GRAPH
#anykey
cap = 1
for cap in range(2):
    print 't={0} E={1}'.format(es.time,es.analysis.energy(es)["ideal"])
#    es.inter(individual = cap) 
    es.non_bonded_inter.set_force_cap(cap)
    es.integrator.run(20)
    i+=19
    

es.non_bonded_inter.set_force_cap(0)

#my_actor = DipolarDirectSumCpu(bjerrum_length = 1.0)
my_actor = DipolarBarnesHutGpu(bjerrum_length = 1.0, epssq = 100.0, itolsq = 4.0)
es.actors.add(my_actor)

for i in range(1000):
    temp = es.analysis.energy()["ideal"] /((deg_free/2.0)*n_part)
    print 't={0} E={1} , T={2}")'.format(es.time,es.analysis.energy()["total"],temp)
    es.integrator.run(10)

#imd positions  


#the option from the Tanygin's manunscript:
#################################inter magnetic 1.0 bh-gpu
#es.inter.add(magnetic, 1.0, bh-gpu)
#inter magnetic 1.0 dds-gpu

#H_demag = 0.0
#dipm_obs = es.observables.MagneticDipoleMoment(com_dipole_moment, all) 

print ("==============================")
print ' {0} | {1} | {2} | {3}'.format('t_epoch','t_in_silico','Mz',T)
print ("==============================")
es.integrator.run(steps = 0,recalc_forces = True) ############?
H = np.array([Hx,Hy,((Hz*Oe_to_Am)/H0)])
visualizer = visualization.mayaviLive(es)

def main_tread():        
    i=0
    timestep = 2
    while True:
        for n in range(n_part):
            es.part[n].ext_torque = np.cross(es.part[n].dip,H)
        temp = es.analysis.energy()["ideal"]/((deg_free/2.0)*n_part)
        #print '(t={0} E={1} , T={2})'.format(es.time,es.analysis.energy()["total"],temp)
      #  print '(t={0})'.format(es.time)
      #  print (str(int(es.time%analyse_step)) == str(int(timestep)))
      #  print (str(int(es.time%analyse_step)))
        if(str(int(es.time%analyse_step)) == str(int(timestep))):        
            print('{0};\n'.format(es.time))
            for k in range(n_part):
                print '{0};{1};{2};{3};\n'.format(k,es.part[k].pos[0],es.part[k].pos[1],es.part[k].pos[2])
        es.integrator.run(steps = 10,recalc_forces = False)
        i=i+10   
        es.integrator.run(10)
        visualizer.update()

t = Thread(target = main_tread)
t.daemon = True
t.start()
visualizer.start()

            
    
#puts "t=[setmd time] E=[analyze energy total], T=$temp"
#puts "t=[setmd time] dipm=[observable $dipm_obs print]"
   # H_demag = 0.0
#set H_demag [expr -[lindex [ observable $dipm_obs print ] 2]/$V_carr ]
   # Mz = [lindex( observable ,dipm_obs ), 2]/V_carr
   #  H_mag = ((H_ext_Oe*Oe_to_Am)/H0+H_demag)
   # es.constraint.add(ext_magn_field, 0, 0, H_mag)
#############################################integrate 10 reuse_forces?
#imd positions 

print '{0},{1},{2},{3}'.format('[clock seconds]',es.time,'Mz',temp)

#imd disconnect 


