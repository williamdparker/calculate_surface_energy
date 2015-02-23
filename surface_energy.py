#!/usr/bin/env python


#Conversion factors
meters_per_bohr = 0.529177 * 10**-10
eV_per_Ry = 27.211385 / 2
J_per_eV = 1.60217657 * 10**-19

#Experimental value of fcc Pt lattice spacing (in Bohr)
lattice_parameter = 7.407720

#Pt 110 surface lattice vectors
lattice_vectors = [[1., 0., 0.],[0., 0.707107, 0.],[0.,0.,7.521092]]

#111 surface lattice vectors
#lattice_vectors = [[0.707107, 0., 0.],[-0.353553, 0.612372, 0.],[0.,0.,8.082904]]

#100 surface lattice vectors
#lattice_vectors = [[0.707107, 0., 0.],[0., 0.707107, 0.],[0.,0.,8.399772]]


##DFT-PBE values
#from literature
# (110)
#1.85  J/m^2 in Singh-Miller and Marzari 80, 235407 (2009)
#1.38 eV/atom (unrelaxed) "   "           "    "  
#1.30 eV/atom   (relaxed) "   "           "    "  
#
# (111)
#1.81  J/m^2 in Singh-Miller and Marzari 80, 235407 (2009)
#0.90 eV/atom (unrelaxed) "   "           "    "  
#
# (100)
#1.49  J/m^2 in Singh-Miller and Marzari 80, 235407 (2009)
#0.65 eV/atom (unrelaxed) "   "           "    "  
#
#my calculations:
#110 slab Pt total energy in Ry per simulation cell (7 atoms)
e_Pt110_surface_unrelaxed = -366.88610122
#fcc Pt total energy in Ry per simulation cell (1 atom)
e_Ptfcc_bulk = -52.44208196
#Number of atoms in slab simulation
natom =  7

##DMC values
#fcc Pt total energy in Ha per simulation cell (1 atom)
e_Ptfcc_bulk_DMC = -25.67956
#111 slab Pt total energy in Ha per simulation cell (16*7 atoms) 
e_Pt111_surface_DMC = -2875.184
de_Pt111_surface_DMC = 0.002
#100 slab Pt total energy in Ha per simulation cell (16*7 atoms)
e_Pt100_surface_DMC = -2874.597
de_Pt100_surface_DMC = 0.025
#Number of atoms in slab simulation
natom_DMC = 16 * 7

def calculate_simulation_cell_surface_area(parameter, vectors):
    area = parameter**2 * (vectors[0][0] * vectors[1][1] - vectors[0][1] * vectors[1][0])
    return area

def calculate_surface_energy(slab_energy,slab_size,bulk_energy_per_atom):
    gamma = 0.5 * (slab_energy - slab_size * bulk_energy_per_atom)
    return gamma

area = calculate_simulation_cell_surface_area(lattice_parameter,lattice_vectors)
area_in_square_meters = area * (meters_per_bohr)**2

print "Simulation cell surface area = {} bohr**2 = {} m**2".format(area,area_in_square_meters)

#Surface energy
gamma = calculate_surface_energy(e_Pt110_surface_unrelaxed,natom,e_Ptfcc_bulk)
gamma_in_eV_per_atom = gamma * eV_per_Ry
gamma_in_J_per_square_meter = gamma_in_eV_per_atom * J_per_eV / area_in_square_meters

print "Gamma = {} Ry/atom = {} eV/atom = {} J/m**2".format(gamma,gamma_in_eV_per_atom,gamma_in_J_per_square_meter)


