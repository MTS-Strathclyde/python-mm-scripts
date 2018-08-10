'''
Created on Apr 4, 2012

@author: maxim
'''
import re
import pybel

class Molecule:
    def __init__(self, filename, low, high):
        '''
        Constructor
        '''
        self.filename = filename
        self.xyz_filename = self.filename + ".xyz"
        self.low = low
        self.high = high
              
    def stripCells(self):
        self.cycle_atom_nums = ()           #tuple with number of atoms with ## in front
        f = open(self.filename, 'r')        #file open
        l = []                              #list with stripped from ## lines
        line_num = 1
        for line in f:
            if "##" in line[:4]:        #due to spaces, sometimes added by geom. generation programs
                self.cycle_atom_nums += (line_num - 3,)     #program controls first 4 elements of line (not 2)
                line = line.strip("# ")     #line is striped from spaces as well
            l.append(line)
            line_num += 1
        #Checking for number of cycles
        if len(self.cycle_atom_nums) != 2:
            raise ValueError("Molecule must contain exactly one cycle!")
        #Creating stripped text version
        self.stripText = "".join(l)
        print 'Checking distance between atom num ' + str(self.cycle_atom_nums[0]) + " and atom num " + str(self.cycle_atom_nums[1])
        f.close()
    
    def edit(self):
        """creates pure z-matrix from given file, ready to be sent to pybel"""
        self.stripCells()
        i = 0
        self.angle0_incr_list = []
        #angle0 incr list contains lists of starting angles and increments
        while True:
            match = re.search(r'(\s\d*\.\d+\.*|\s\d+\.*)!(\d*\.\d+\.*|\d+\.*\s)', self.stripText)
            #r - ignore \
            #\s must be space befor num!num and after it
            #\d - digit
            # * - 0 or more times
            # + - 1 or more times
            # d*\.\d+.*  - float, which may start as .343 and end with .
            # \d+.* - integer, which may end with .
            
            if match != None:
                l = [float((match.group(1)).strip(" .")), float((match.group(2)).strip(" ."))]
                self.angle0_incr_list.append(l)
                self.stripText = re.sub(r'(\s\d*\.\d+\.*|\s\d+\.*)!(\d*\.\d+\.*|\d+\.*\s)', ' {0[' + str(i) + ']} ', self.stripText, 1)
                i += 1
            else:
                if i == 0:
                    print("0 conformers will be created from this file!")
                self.text_for_formating = self.stripText
                print "List of start angles and increments:", self.angle0_incr_list
                break

    def zMatrix(self, angle_list):
        """Creates zMatrix from given state and returns it as string"""        
        return self.text_for_formating.format(angle_list)
            
    def goalTest(self, angle_list, state_num):
        zMatrix = self.zMatrix(angle_list)
        pybel_mol = pybel.readstring("mopin", zMatrix)
        atoms = pybel_mol.atoms
        atom1 = atoms[self.cycle_atom_nums[0] - 1].OBAtom
        atom2 = atoms[self.cycle_atom_nums[1] - 1].OBAtom
        d = atom1.GetDistance(atom2)
        accepted = self.high > d > self.low
        if accepted:
            #2 possibility for writing output:
            #first with pybel
            #it is faster, but babel modifies z - matrix after transform, which might be bad
            #also he removes keywords
            
            #second using python
            #it is slower, but doesn't change z-mat and removes keywords
            
            #pybel_mol.write(format='mopin', filename = str(state_num) + '_' + self.filename)
            f = open(str(state_num) + '_' + self.filename, 'w')
            f.write(zMatrix)
            f.close()
        return accepted
        

if __name__ == '__main__':
    ##Constructor test
#    m = Molecule("bata")
#    print m.filename, m.xyz_filename
    ##stripCells test
#    m = Molecule("01_zmat_test.dat")
#    m.stripCells()
#    print(m.stripText)
#    print
    ##edit test
    m = Molecule("01_zmat_test.dat", 1, 3)
    m.edit()
    print(m.text_for_formating)
    ##zmat test
    print m.zMatrix("111")
    ##openbabel test
    mol = pybel.readstring("mopin", m.zMatrix("111"))
    atoms = mol.atoms
    atom1 = atoms[m.cycle_atom_nums[0]].OBAtom
    atom2 = atoms[m.cycle_atom_nums[1]].OBAtom
    print atom1.GetDistance(atom2)
    print atom2.GetDistance(atom1)
    print atom1.GetDistance(atom1)
    print m.goalTest("111")
    
    
    
    
    
    
    
    
    