#!/usr/bin/env python2
import sys
import numpy as np


K_B = 8.314/4184

RESULTS = """dGsolv(closure)= {exchem} kcal/mol
PMV= {pmv} A^3

dGsolv(PC+)= {PC_plus} kcal/mol
dGsolv(PC)= {PC} kcal/mol

P_minus_ideal_gas_pressure= {pressure_plus} kcal/mol/A^3
P= {pressure} kcal/mol/A^3

"""



class Xvv(object):
    """ Wrapper around xvvfile used to compute 3d-rism pressure """
    def __init__(self, fname):
        """ Read xvvfile and set instance attributes 
        
        Parameters
        ----------
        
        fname : string
            Path to a valid xvv file
        """
        self.fname = fname
        self.ngrid = None
        self.nsites = None
        self.nspecies = None
        self.temperature = None
        self.dr = None
        self.atom_names = None
        self.densities = None
        self.xvv_data = None
        self.multiplicities = None
        self.unique_sites_per_species = None
        self.total_sites_per_species = None
        self.compressibility = None
        self.species_densities = None
        self.normalized_densities = None
        self._read_xvvfile()
        self._compute_species_properties()

    def _read_xvvfile(self):
        with open(self.fname) as f:
            lines = f.readlines()
        tot_lines = len(lines)
        for i, line in enumerate(lines):
            line = line.split()
            if len(line) <= 1:
                continue
            if line[1] == 'POINTERS':
                data = map(int, lines[i+2].split())
                self.ngrid, self.nsites, self.nspecies = data
            if line[1] == 'MTV':
                self.multiplicities = map(int, lines[i+2].split())
            if line[1] == 'NVSP':
                self.unique_sites_per_species = map(int, lines[i+2].split())
            if line[1] == 'THERMO':
                data = lines[i+2].split()
                self.temperature = float(data[0]) # K
                self.dr = float(data[4]) # Angstrom
                self.compressibility = float(data[3]) # Angstrom^{3}
            if line[1] == 'ATOM_NAME':
                data = lines[i+2].strip()
                #split into groups of 4
                self.atom_names = [data[i:i+4].strip() for i in range(0, len(data), 4)]
            if line[1] == 'RHOV' and len(line) == 2:
                self.densities = map(float, lines[i+2].split())
                #are there more lines with density?
                counter = 3
                while lines[i+counter].startswith(' '):
                    self.densities.extend(map(float, lines[i+counter].split()))
                    counter += 1
                try:
                    assert len(self.densities) == len(self.atom_names)
                except AssertionError:
                    print('Inconsistent number of densities and atom names')
                    print(self.densities)
                    print(self.atom_names)
                    raise ValueError
            if line[1] == 'XVV' and len(line) == 2:
                self.xvv_data = []
                xvv_ind = i + 2
                while xvv_ind < tot_lines and not lines[xvv_ind].startswith('%'):
                    self.xvv_data.extend(lines[xvv_ind].split())
                    xvv_ind += 1
                break
        assert len(self.xvv_data) == self.ngrid*self.nsites*self.nsites
        self.xvv_data = np.array(self.xvv_data, dtype=float)
        self.xvv_data = np.reshape(self.xvv_data,
                                   (self.ngrid, self.nsites, self.nsites),
                                   order='F')

    def _compute_species_properties(self):
        self.normalized_densities = []
        for density, multiplicity in zip(self.densities, self.multiplicities):
            self.normalized_densities.append(density/multiplicity)
        self.species_densities = []
        self.total_sites_per_species = []
        pointer = 0 
        for sp_sites in self.unique_sites_per_species:
            pointer += sp_sites
            total_sites = sum(self.multiplicities[pointer - sp_sites:pointer])
            self.total_sites_per_species.append(total_sites)
            self.species_densities.append(self.normalized_densities[pointer - 1])
        assert len(self.species_densities) == self.nspecies
    
    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. The pressure can be pretty
            sensitive to it. It is recommended to experiment with a couple of
            k values or better, plot dependency of pressure on it to see
            which value works best.
            
        Return
        ------
        pressures : tuple of floats
            Tuple containeing two pressures.
            First element is 3D-RISM pressure (used in PC), second element is
            3D-RISM pressure minus ideal gas pressure (used in PC+).
            Both have units of kcal/mol/A^3.
        """
        xvv_k = self.xvv_data[k,:,:]
        density_vec = np.array(self.normalized_densities)
        mult_vec = np.array(self.multiplicities)
        # Z_k from sergievskyi's article
        z_k = mult_vec/density_vec*(np.identity(self.nsites) - np.linalg.inv(xvv_k))
        z_k_sum_densities2 = np.sum(density_vec*z_k*density_vec.T)
        densities_times_sites = [sites*dens for sites, dens in zip(self.total_sites_per_species,
                                                                   self.species_densities)]
        pressure = sum(densities_times_sites) - .5*z_k_sum_densities2
        pressure = pressure*self.temperature*K_B
        ideal_pressure  = sum(self.species_densities)*K_B*self.temperature
        return pressure, pressure - ideal_pressure



RESULTS_NAME = 'results2.txt'

def main(log):
    with open(log) as f:
        txt = f.read()
    for l in txt.splitlines():
        if 'xvvfile=' in l:
            xvv_p = l.split('=')[1]
            xvv = Xvv(xvv_p)
        if l.startswith('rism_excessChemicalPotential'):
            ex = float(l.split()[1])
        if l.startswith('rism_partialMolarVolume'):
            pmv = float(l.split()[1])
        if l.startswith('rism_KirkwoodBuff'):
            kbs = l.split()[1:]
            kbs = np.array(map(float, kbs))
        if l[0:11] == "rism_exchem":
            ex = float(l.split()[1])
        if l[0:11] == "rism_volume":
            pmv = float(l.split()[1])
        if l.startswith('rism_KB'):
            kbs = l.split()[1:]
            kbs = np.array(map(float, kbs))
    vol_ex = pmv - xvv.compressibility
    p, p_plus = xvv.compute_3drism_pressures()
    pc = ex - vol_ex*p
    pcp = ex - vol_ex*p_plus
    #print pc, pcp
    results = RESULTS.format(exchem=ex, pmv=vol_ex, PC=pc, PC_plus=pcp,
                             pressure=p, pressure_plus=p_plus)
    with open(RESULTS_NAME, 'w') as f:
        f.write(results)
        
    
if __name__ =='__main__':
    main(sys.argv[1])
    
    
            
            
            
            
            
            


