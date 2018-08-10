import csv
import sys

"""
To start script:
python Separate_by_energies.py [increment_bound] file.csv

Takes as input csv file with sorted energies and molecule names in form:
mol_name, energy
mol_name2, energy2
...
where energy2 > energy

Divides molecules into groups with simillar energy. Group can be arbitrarily
big. As long as next molecule energy is smaller then incr_bound (default is 0.1)"""

def divide_by_en(incr_bound, csv_f):
    """(incr_bound, csv_f)
    Produces list of groups, which contain filenames of molecules divided by energies    
    """
    #Make file object from suplied filename
    f = file(csv_f)
    #Initialize csv reader
    reader = csv.reader(f)
    #Create list of groups
    groups = []
    #Previous energy
    pr_en = -float("inf")
    #Main cycle
    for line in reader:
        #energy of this molecule
        ener = float(line[1])
        #If this energy is bigger, then maximum allowed increment and 
        #previous group is bigger then 1
        if ener - pr_en > incr_bound:
            #Create new group
            groups.append([])
            #And add there this molecule name
            groups[-1].append(line[0])
            #Print info
            if len(groups) > 1:
                print "Group ends with energy:",pr_en
                print "It has",len(groups[-2]),"molecules."
            print "New group starts with energy:",ener
        #Else
        else:
            #Simply append current molecule to existing group
            groups[-1].append(line[0])
        #Update previous energy
        pr_en = ener
    #remove groups with one molecule
    corrected_groups = [groups[0]]
    for i in range(1, len(groups)):
        if len(groups[i]) < 2:
            corrected_groups[-1] += groups[i]
        else:
            corrected_groups.append(groups[i])
    #Info about last group
    print "Group ends with energy:",pr_en
    #Close file
    f.close()
    return corrected_groups

def is_float(num):
    """Checks, weather argument is float. Accepts only strings and numbers"""
    try:
        float(num)
        return True
    except ValueError:
        return False 

def main(argv):
    """
    [[-increment], CSV_file]
    Separates by energies. Algorithm specified in script description.
    """
    #Default increment is 0.1
    incr_bound = 0.1
    #If user suplied another as first argument in lunching program, it should
    #be changed
    if is_float(argv[0]):
        incr_bound = float(argv[0])
        csv_f = argv[1]
    else:   #First argument is csv file
        csv_f = argv[0]
    #Now call the function; which separates by energies.
    #It returns list of different list, each containing filnames of molecules
    #belonging to the same group.
    return divide_by_en(incr_bound, csv_f)

if __name__=="__main__":
    main(sys.argv[1:])
