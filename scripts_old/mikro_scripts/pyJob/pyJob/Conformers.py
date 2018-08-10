'''
Created on Apr 25, 2012

@author: maxim
'''

from math import ceil


def bounds_gen(angle0_incr_list):
    """how many multiples of given angles will be there"""
    bounds = []
    for index in range(len(angle0_incr_list)):
        angle0 = float(angle0_incr_list[index][0])
        incr = float(angle0_incr_list[index][1])
        nurgad_kokku = ceil((360 - angle0)/(incr))
        bounds.append(int(nurgad_kokku))
    return bounds


def conformers(molecula):
    angle0_incr_list = molecula.angle0_incr_list
    bounds = bounds_gen(angle0_incr_list)
    total_angles = len(angle0_incr_list)

    print("angle0_incr_list", angle0_incr_list)
    print("bounds", bounds)

    #global var for numeration
    global state_num
    state_num = 1
    
    accepted = []
    
    def abi(angle_list, index):
        """recursive function for angle list generation"""
        if index < total_angles:
            angle0 = angle0_incr_list[index][0]
            incr = angle0_incr_list[index][1]
            for multiple in range(bounds[index]):
                copy = angle_list[:]
                copy.append(angle0 + incr*multiple)
                
                if len(copy) == total_angles:
                    global state_num
                    #print(copy, state_num) for test
                    if molecula.goalTest(copy, state_num):
                        accepted.append(str(state_num) + '_' + molecula.filename)
                    state_num += 1
                abi(copy, index + 1)
    abi([], 0)
    return accepted

if __name__ == '__main__':
    class test_class:
        def __init__(self, angle0_incr_list):
            self.angle0_incr_list = angle0_incr_list
    
    molec = test_class([[0, 60], [0, 60], [0, 60], [0, 60]])
    conformers(molec)
