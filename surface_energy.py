#!/usr/bin/env python


#Conversion factors
meters_per_bohr = 0.529177 * 10**-10
eV_per_Ry = 27.211385 / 2
J_per_eV = 1.60217657 * 10**-19

#Experimental value of fcc Pt lattice spacing (in Bohr)
lattice_parameter = 7.407720

lattice_vectors = {
        110:[[1., 0., 0.],[0., 0.707107, 0.],[0.,0.,7.521092]],
        111:[[0.707107, 0., 0.],[-0.353553, 0.612372, 0.],[0.,0.,8.082904]],
        100:[[0.707107, 0., 0.],[0., 0.707107, 0.],[0.,0.,8.399772]]
        }

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
total_energies = {
        'dft':{
            'fcc':[-52.44208196, 0, 1],
            '110':[-366.88610122, 0, 7]
            },
        'qmc':{
            'fcc':[-25.67956, 0, 1],
            '111':[-2875.184, 0.002, 112],
            '100':[-2874.597, 0.025, 112]
            }
        }

def calculate_simulation_cell_surface_area(parameter, vectors):
    area = parameter**2 * (vectors[0][0] * vectors[1][1] - vectors[0][1] * vectors[1][0])
    return area

def calculate_surface_energy(slab_energy,slab_size,bulk_energy_per_atom):
    gamma = 0.5 * (slab_energy - slab_size * bulk_energy_per_atom)
    return gamma

area = calculate_simulation_cell_surface_area(lattice_parameter,lattice_vectors[110])
area_in_square_meters = area * (meters_per_bohr)**2

print "Simulation cell surface area = {} bohr**2 = {} m**2".format(area,area_in_square_meters)

#Surface energy
slab_energy = total_energies['dft']['110'][0]
slab_size = total_energies['dft']['110'][2]
bulk_energy = total_energies['dft']['fcc'][0]
gamma = calculate_surface_energy(slab_energy, slab_size, bulk_energy)
gamma_in_eV_per_atom = gamma * eV_per_Ry
gamma_in_J_per_square_meter = gamma_in_eV_per_atom * J_per_eV / area_in_square_meters

print "Gamma = {} Ry/atom = {} eV/atom = {} J/m**2".format(gamma,gamma_in_eV_per_atom,gamma_in_J_per_square_meter)


