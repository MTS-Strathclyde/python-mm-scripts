#!/usr/bin/python2.7
import sys
import os
import random
import shutil

#divides randomly files in given folder into n groups and puts them in separate
#folders

def copy_to_folder(mol_list, i):
    """copies given molecule names to folder named _present_folder_name_sub_i"""
    #Get folder name
    folder_name = os.getcwd().split("/")[-1] + "_sub_" + str(i)
    #Make folder
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    else:
        folder_name += "c"
        os.makedirs(folder_name)
    #Copy
    for filename in mol_list:
        shutil.copy(filename, os.getcwd() + "/" + folder_name)    

def main():
    if sys.argv[1].isdigit():
        n = int(sys.argv[1])
        files = sys.argv[2:]
    else:
        n = 2 #by default devides into two folders
        files = sys.argv[1:]
   ## print "files", files        
    random.shuffle(files)   #shuffle files
   ## print "shuffled", files
    num_f_sub = len(files)/n  #amount of files in subset
    for i in range(n - 1):  #loop
        sub_set = files[i*num_f_sub:(i+1)*num_f_sub]
      ##  print "copying this:", sub_set
        copy_to_folder(sub_set, i)
    sub_set = files[(n-1)*num_f_sub:]
   ## print "copying this:", sub_set
    copy_to_folder(sub_set, n-1)
    
if __name__ == "__main__":
    main()
