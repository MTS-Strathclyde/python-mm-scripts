import sys

"""Reads given arc file and extracts optimized torison angles from it"""

def get_torisons(arc_f):
    #convert name to file object
    f = file(arc_f)
    #line number
    l = 0
    #list of dihedrals
    dihed = []
    #loop over all lines in file
    for line in f:
        l += 1 #increment line number
        if l > 45: #Z-mat starts from line 46
            #split line into list of strings
            strings = line.split()
            #if this is not last line
            if len(strings) > 1:
                #dihedral is 5-th element in produced list
                dihedral = float(strings[5])
                #add to list
                dihed.append(dihedral)
    f.close()
    return dihed

def main():
    #.arc file is passed as argument
    arc_f = sys.argv[1]
    #Function returns list of dihedral angles
    print get_torisons(arc_f)

if __name__ == "__main__":
    main()
