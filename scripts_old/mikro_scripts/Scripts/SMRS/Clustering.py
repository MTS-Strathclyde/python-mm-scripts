import pybel
import openbabel as ob
import sys
import os
import csv

"""

Evaluates representative molecules
To start:
python Clustering.py rep1.out [rep2.out rep3.out ...]
(in folder, which contatains molecules to be clustered)

RESULT:
Prints number of molecules, which belong to different representative molecules and average
rmsd between them

TODO:
This script should save CSV with ditailed statistics. Probably matrix of all representatives and RMSD's
Also there should be more clear way to avaluate fittness of representatives.
Better to be one number score.

"""

def cluster(ob_rep_mols):
    """Function creates dictionary, which keys are representive molecules and values are
dictionaries, which contain molecule names as keys and RMSD scores as values."""
    #Firstly dictionary, which represents clusters is created.
    #It's done by zipping 2 lists - one with rep. molecules, other with empty dictionaries
    dic = dict(zip(ob_rep_mols, [{} for i in ob_rep_mols]))
    #Then OBAlign object is initialized
    align = ob.OBAlign(False, False)
    #After that function loops over all files in directory
    for f in os.listdir("."):
        #If file name ends with .out, it indicates, that it is mopout file
        if f[-4:] == ".out":
            #Function converts it and uses as reference in OBAlign object
            align.SetRefMol(pybel.readfile("mopout", f).next().OBMol)
            #Then we need to find, to which cluster belongs given molecule
            #For that we setup two vaiables
            ref_m = None    #Ref molecule, to which molecule in dir belongs
            rmsd_m = 1000     #And RMSD score between them
            #And start loop
            for ref in ob_rep_mols:
                #It loops over representative molecules and set ups them as target molecules
                #in OBAling file
                align.SetTargetMol(ref)
                #Then it checks RMSD between them (if it is possible)
                if align.Align():
                    rmsd = align.GetRMSD()
                    if rmsd < rmsd_m:
                        rmsd_m = rmsd
                        ref_m = ref
            #Inner loop end
            #After end of loop molecule in directory goes to dictionary as the value of that
            #representative molecule, to which it is more similair
            #RMSD between them is also saved
            print "Min rmsd:", rmsd_m
            dic[ref_m][f] = rmsd_m
    #Loop over all files in directory ends
    #Dictionary is populated and returned
    return dic

def main():
    #Representative molecule names are entered as arguments
    rep_mols = sys.argv[1:]
    #Later they are converted into list of OBMolecules
    ob_rep_mols = [pybel.readfile("mopout", mol).next().OBMol for mol in rep_mols]
    #And passed to clustering function, which returns dictionary
    dic = cluster(ob_rep_mols)
    #For now I will simply print average RMSD's and number of mols in clusters
    for rep in dic.keys():
        print "Representative molecule:",rep.GetTitle()
        print "Has",len(dic[rep]),"conformers belonging to it."
        print "And average RMSD:",sum(dic[rep].values())/len(dic[rep])
        
if __name__ == "__main__":
    main()
