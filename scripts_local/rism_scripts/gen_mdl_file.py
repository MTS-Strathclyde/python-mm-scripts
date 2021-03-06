#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 17:08:12 2015

@author: max
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import os
import datetime
import sys
import argparse
import subprocess
import glob
import distutils.spawn



PREPDAT_PL_SCRIPT_BEGIN = """#!/usr/bin/perl -w
use POSIX;
use strict;

#This script creates parameter files to be read by RISM1D for all predefined solvent 
#and ion parameters in the Amber force field (FF99SB is used).
#This is done in two parts for each molecule.  First, a PRMTOP and CRD file are 
#created in TLEAP. These two resulting files are parsed and the relavent
#data placed into a single file following the PRMTOP format and can be read 
#by NXTSEC.

#detect atoms with the same sigma, epsilon and charge and assign the correct 
#multiplicity value.  This makes for a small number of solvent atoms and makes 
#3D-RISM more efficient
my $multiplicity=1;

#use the Amber PARM7 format.  This is really the only supported format in RISM1D.
#The other format is namelist but this is difficult to read for variable sized arrays
my $ambParm=1;
"""
PREPDAT_PL_SCRIPT_CUSTOM_PART = """#solvent models
my %solvent=(
	  "{name}"=>["{solvent_type}","{name}"]);
"""
PREPDAT_PL_SCRIPT_END ="""#output directory
my $datdir=".";

#proceed through each solvent
foreach my $solvent (keys (%solvent)) {
    print "solvent $solvent\n";
    tleap($solvent,@{$solvent{$solvent}},$datdir);
    buildParm($solvent,$datdir);
}

exit;

#creates a single residue system in TLEAP for the given molecule.
#$solvent - residue name
#$type    - i.e. WAT, ION, OTH...
#$alt     - depricated in this version, same as name
#           prevously: (frcmod suffix (typically not needed for IONs))
#$datdir  - destination for resulting PARM7 and RST7 files
sub tleap {
    my ($solvent,$type,$alt,$datdir) = @_;

    my $tleap = "tleap";
    my $inp = "$datdir/$solvent.inp";
    my $parm = "$datdir/$solvent.parm7";
    my $rst = "$datdir/$solvent.rst7";
    my $log = "$datdir/$solvent.log";
    open(INP,">$inp") || die "ERROR: could not open >$inp:$!";
    print INP "source leaprc.gaff";
    if($type eq "WAT"){
	print INP "WAT = $solvent";
	if($solvent ne "TP3"){
	    print INP "loadAmberParams ".(lc($alt)).".frcmod";
	}
    }elsif($type eq "OTH"){
	    print INP "loadAmberParams ".(lc($alt)).".frcmod\n";
	    print INP "loadAmberPrep ".lc($alt).".prepin";
    }elsif($type eq "ION"){
	
    }else{
	die "ERROR: unknown solvent type: $solvent $type";
    }
    #print INP "sol = sequence { $solvent }";
    print INP "saveAmberParm MOL $parm $rst\n";
    print INP "quit";
    close INP;

    my $cmd="$tleap -f $inp >$log 2>&1";
    !system($cmd) || die "ERROR: could not run $cmd:$!";
}


#create the MDL file for the molecule.  A PARM7 and RST7 file should exist in 
#$datdir.  These are read in and processed for output in the MDL file and then 
#deleted.
#$solvent - solvent molecule name.  Used as the root name for all files
#$datdir  - output location for all work
sub buildParm {
     my ($solvent,$datdir) = @_;
     local $_;
     my $rst="$solvent.rst7";
     my $parm="$solvent.parm7";
     my $mdl="$solvent.mdl";
     #all the parameters we will write out
     my (@coord, @name, @charge, @mass, @typeNum, @typeName, @lja, @ljb, 
	 @radius, @ljIndex, @ljsigma, @ljepsilon,@multi);
     my $i=0;
     #get coordinates
     open(RST,"<$rst") || die "ERROR: could not open <$rst:$!\n";
     <RST>;
     <RST>;
     while(<RST>){
	 
	 my @temp=split();
	 push(@coord, [@temp[0,1,2]]);	 
	 $i++;
	 if($#temp==5){
	     push(@coord, [@temp[3,4,5]]);
	     $i++;	     
	 }
     }
     close RST;

     #read all parameters into a hash
     my %top = readPRMTOP($parm);

     #compute the well depth and radius for each atom.
     my $nType = sqrt($#{$top{NONBONDED_PARM_INDEX}}+1);
     for(my $i=0;$i<=$#{$top{ATOM_TYPE_INDEX}};$i++){
	 my $index = ${$top{NONBONDED_PARM_INDEX}}[$nType*(${$top{ATOM_TYPE_INDEX}}[$i]-1) + ${$top{ATOM_TYPE_INDEX}}[$i]-1]-1;
	 if($index < 0 || ${$top{LENNARD_JONES_BCOEF}}[$index] == 0 || ${$top{LENNARD_JONES_ACOEF}}[$index] == 0){
#	     push(@ljepsilon, 0.046);
#	     push(@ljsigma, .2*pow(2,1./6.));
	     #if the radius is zero, assign a sero for further processing
	     push(@ljepsilon, 0);
	     push(@ljsigma, 0);
	     next;
	 }
	 push(@ljepsilon, ${$top{LENNARD_JONES_BCOEF}}[$index]*${$top{LENNARD_JONES_BCOEF}}[$index]/(4*${$top{LENNARD_JONES_ACOEF}}[$index]));
	 push(@ljsigma, pow(2*${$top{LENNARD_JONES_ACOEF}}[$index]/(${$top{LENNARD_JONES_BCOEF}}[$index]),1./6.)/2.);
     }

     #use the coincident radius for atoms that have none defined
     for(my $index = 0; $index <= $#ljsigma; $index++){
	 !$ljsigma[$index] || next;
	 my @parent = getParents(\%top,$index);
	 my @bondlength;
	 for my $parent (@parent){
	     push(@bondlength,getBondLength(\%top,$index,$parent));
	 }
	 ($ljsigma[$index],$ljepsilon[$index]) = calcCoincidentRadius(\@parent,\@bondlength,\@ljsigma,\@ljepsilon);
	 #Radius must be greater than 0
	 $ljsigma[$index] >0 || warn "Negative coincident radius for atom ".($index+1).
	     "${$top{ATOM_NAME}}[$index]\n";
     }
     
     if($multiplicity){
	 removeMultiples(\@coord, \@{$top{ATOM_NAME}}, \@{$top{CHARGE}}, \@{$top{MASS}}, \@{$top{ATOM_TYPE_INDEX}}, \@{$top{AMBER_ATOM_TYPE}}, \@{$top{RADII}}, \@{$top{NONBONDED_PARM_INDEX}}, \@ljsigma, \@ljepsilon, \@multi);
     }

     #write the output in PARM7 or namelist format
     open(MDL,">$solvent.mdl") || die "ERROR: could not open >$solvent.mdl:$!\n";
     if($ambParm){
	 writeAmbHeader(*MDL);
	 writeAmbParm(*MDL,"20a4","TITLE",[$solvent]);
     }else{
	 print MDL "&SOLVENT\n";
     }

     #write array sizes
     if($ambParm){
	 my @sizes = ($#coord+1,$#{$top{ATOM_TYPE_INDEX}}+1);
	 writeAmbParm(*MDL,"10I8","POINTERS",\@sizes);
     }

     #write atom types
     if($ambParm){
	 writeAmbParm(*MDL,"10I8","ATMTYP",\@{$top{ATOM_TYPE_INDEX}});
     }else{
	 writeParam(*MDL,"!ATOM TYPES","ATMTYP",\@{$top{ATOM_TYPE_INDEX}});
     }

     #write atom names
     if($ambParm){
	 writeAmbParm(*MDL,"20a4","ATMNAME",\@{$top{ATOM_NAME}});
     }else{
	 writeParam(*MDL,"!ATOM NAMES","ATMNAME",\@{$top{ATOM_NAME}});
     }

     #write masses
     if($ambParm){
	 writeAmbParm(*MDL,"5e16.8","MASS",\@{$top{MASS}});
     }else{
	 writeParam(*MDL,"!MASS","MASS",\@{$top{MASS}});
     }

     #write charges
     if($ambParm){
	 writeAmbParm(*MDL,"5e16.8","CHG",\@{$top{CHARGE}});
     }else{
	 writeParam(*MDL,"!CHARGES","CHG",\@{$top{CHARGE}});
     }

     #write LJ epsilon parameters
     if($ambParm){
	 writeAmbParm(*MDL,"5e16.8","LJEPSILON",\@ljepsilon);
     }else{
	 writeParam(*MDL,"!LENNARD-JONES EPSILON PARAMETER","LJEPSILON",\@ljepsilon);
     }

     #write LJ sigma parameters
     if($ambParm){
	 writeAmbParm(*MDL,"5e16.8","LJSIGMA",\@ljsigma);
     }else{
	 writeParam(*MDL,"!LENNARD-JONES SIGMA PARAMETER","LJSIGMA",\@ljsigma);
     }
     #write multiplicity
     if($multiplicity){
	 if($ambParm){
	     writeAmbParm(*MDL,"10I8","MULTI",\@multi);
	 }else{
	     writeParam(*MDL,"!ATOM TYPE MULTIPLICITY","MULTI",\@multi);
	 }
     }
     
     
     #write coordinates
     if($ambParm){
	 my @flatcoord;
	 for(my $i=0; $i<=$#coord; $i++){
	     for(my $j=0; $j<=$#{$coord[$i]}; $j++){
		 push(@flatcoord,${$coord[$i]}[$j]);
	     }
	 }
	 writeAmbParm(*MDL,"5e16.8","COORD",\@flatcoord);
     }else{
	 print MDL "/\n";
	 for(my $i=0;$i<=$#coord;$i++){
	     print MDL "@{$coord[$i]}\n";
	 }
     }
     close MDL;
     unlink "$solvent.parm7";
     unlink "$solvent.rst7";
}

#looks for atoms with identical charge and Lennard-Jones parameters.  
#If found, the later is deleted and the multiplicity of the remainder is incremented.
sub removeMultiples {
    my ($coord, $name, $charge, $mass, $typeNum, $typeName, $radius, $ljIndex, $ljsigma, $ljepsilon, $multi) = @_;
    for(my $i=0; $i<=$#{$name}; $i++){
	@$multi[$i]=1;
    }
    for(my $i=0; $i<=$#{$name}; $i++){
	for(my $j=$i+1; $j<=$#{$name}; $j++){
	    if(@$ljsigma[$i] == @$ljsigma[$j] && @$ljepsilon[$i] == @$ljepsilon[$j] &&
	       @$charge[$i] == @$charge[$j]){
		my $coordi=0;
		for(my $k=0;$k<=$i;$k++){
		    $coordi+=@$multi[$k];
		}
		my $coordj=0;
		for(my $k=0;$k<=$j;$k++){
		    $coordj+=@$multi[$k];
		}
		$coordj--;
		(@$coord[$coordi],@$coord[$coordj])=(@$coord[$coordj],@$coord[$coordi]);
		splice(@$name,$j,1);
		splice(@$charge,$j,1);
		splice(@$mass,$j,1);
		splice(@$typeNum,$j,1);
		splice(@$typeName,$j,1);
		splice(@$radius,$j,1);
		splice(@$ljIndex,$j,1);
		splice(@$ljsigma,$j,1);
		splice(@$ljepsilon,$j,1);
		splice(@$multi,$j,1);
		@$multi[$i]++;
		$j--;
	    }
	}
    }
}

#reads in the given PARM7 file into a hash using the %FLAG for each section as 
#the key.  Each key then contains a 1D array of values.  The entire hash is returned.
sub readPRMTOP {
  my $prmtop = shift;
  my %prmHash;
  open(PRMTOP,"<$prmtop") || die "ERROR: could not open $prmtop:$!\n";
  my ($flag,$format,$size);
  while(<PRMTOP>){
    if(/^\%VERSION/){
      
    }elsif(/^\%FLAG/){
      chomp;
      m/\%FLAG ([^\s]+)/;
      $flag = $1;
    }elsif(/^\%FORMAT/){
      chomp;
      m/\%FORMAT\(([^\s]+)\)/;
      $format = $1;
      $format =~ m/[A-Za-z](\d+)/;
      $size = $1;
    }else{
      my @data = split(/(.{$size})/);
      for(my $idata=1;$idata <=$#data;$idata+=2){
	push(@{$prmHash{$flag}},$data[$idata]);
      }
    }
  }
  close PRMTOP;
  return %prmHash;
}

#writes a parameter section using the namelist format
#$mdl     - filehandle
#$comment - comment to preceed
#$varName - name of variable to be written
#$data    - 1D array of data values
sub writeParam {
    my ($mdl,$comment,$varName,$data)=@_;
    print $mdl "$comment\n";
    print $mdl "$varName = ";
    foreach my $data (@$data){
	print $mdl "$data, ";
    }
    print $mdl "\n";
}

#writes the PARM7 format header
#$mdl - filehandle
sub writeAmbHeader {
    my ($mdl) = @_;
    my ($sec, $min, $hr, $day, $mon, $year) = localtime;
    $mon++;
    $year = sprintf("%02d", $year % 100);
    print $mdl "%";
    printf $mdl "VERSION  VERSION_STAMP = V0001.000  DATE = %02d/%02d/%02d  %02d:%02d:%02d\n"
	,$mon,$day,$year,$hr,$min,$sec;
}

#writes a parameter section for the PARM7 format
#$mdl     - filehandle
#$format  - string conversion format for reading/writing
#$varname - name of variable to be written
#$data    - 1D array of data to be written
sub writeAmbParm {
    my ($mdl,$format,$varName,$data)=@_;
    print $mdl "\%FLAG $varName\n";
    print $mdl "\%FORMAT($format)\n";
    $format =~ m/(\d+)([A-Za-z]+)([\d\.]+)/;
    my ($col,$type,$precision)=($1,$2,$3);
    $type = lc($type);
    $type=~s/a/s/;
    for(my $i=0; $i<=$#{$data}; $i+=$col){
	for(my $j=0; $i+$j<=$#{$data} && $j < $col; $j++){
	    printf $mdl "%$precision$type",@$data[$i+$j];
	}
	print $mdl "\n";
    }

}

#gets 'parents' or other atoms bonded to $index.  return as an array of indices 
#with 0 offset
#$top   - PARM7 hash
#$index - atom index to get the bonded neighbours of
sub getParents {
    my $top = shift;
    my $index = shift;
    my @parent;
    #find all of the entrys in BONDS_WITHOUT_HYDROGEN and
    #BONDS_INC_HYDROGEN that contain this index
    for(my $ibond=0; $ibond <= $#{$top->{BONDS_INC_HYDROGEN}}; $ibond+=3){
	if(${$top->{BONDS_INC_HYDROGEN}}[$ibond] == $index*3){
	    push(@parent,${$top->{BONDS_INC_HYDROGEN}}[$ibond+1]/3);
	}
	if(${$top->{BONDS_INC_HYDROGEN}}[$ibond+1] == $index*3){
	    push(@parent,${$top->{BONDS_INC_HYDROGEN}}[$ibond]/3);
	}
    }
    for(my $ibond=0; $ibond <= $#{$top->{BONDS_WITHOUT_HYDROGEN}}; $ibond+=3){
	if(${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond] == $index*3){
	    push(@parent,${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond+1]/3);
	}
	if(${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond+1] == $index*3){
	    push(@parent,${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond]/3);
	}
    }
    return @parent;
}

#given two bonded atoms, return the equilibrium bond length
#$top - PARM7 hash
#$A   - atom 1 index
#$B   - atom 2 index
sub getBondLength {
    my $top = shift;
    my $A = shift;
    my $B = shift;
    for(my $ibond=0; $ibond <= $#{$top->{BONDS_INC_HYDROGEN}}; $ibond+=3){
	if((${$top->{BONDS_INC_HYDROGEN}}[$ibond] == $A*3 && 
	   ${$top->{BONDS_INC_HYDROGEN}}[$ibond+1] == $B*3) ||
	   (${$top->{BONDS_INC_HYDROGEN}}[$ibond] == $B*3 && 
	   ${$top->{BONDS_INC_HYDROGEN}}[$ibond+1] == $A*3)){
	    return ${$top->{BOND_EQUIL_VALUE}}[${$top->{BONDS_INC_HYDROGEN}}[$ibond+2]-1];
	}
    }
    for(my $ibond=0; $ibond <= $#{$top->{BONDS_WITHOUT_HYDROGEN}}; $ibond+=3){
	if((${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond] == $A*3 && 
	   ${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond+1] == $B*3) ||
	   (${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond] == $B*3 && 
	   ${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond+1] == $A*3)){
	    return ${$top->{BOND_EQUIL_VALUE}}[${$top->{BONDS_WITHOUT_HYDROGEN}}[$ibond+2]-1];
	}
    }
    #atoms apparently are not bonded
    warn "bond not found $A $B\n";
}

#given a number of atoms and there bondlengths to the atom in questions, calculates
#the coincident radius.  This is the radius such that LJ potential of our enclosed
#atom is 0 at the same radial distance from the enclosing atom along the bond.  
#I.e. the enclosing radius - the bond length.  This can be a negative number if the
#target atom is not actually enclosed.  This sub routine picks the largest coincident
#radius and returns it along with a epsilon value one tenth that of the enclosing
#atoms.
#$parents     - indices of atoms bonded to the target atom
#$bondlengths - equilibrium bond lengths between parent atoms and the target
#$ljsigma     - list of all sigma values so far determined
#$ljepsilon   - list of all epsilon values so far determined
sub calcCoincidentRadius {
    my $parents = shift;
    my $bondlengths = shift;
    my $ljsigma = shift;
    my $ljepsilon = shift;
    
    my @radius;
    for(my $i = 0; $i<= $#{$parents}; $i++){
	push(@radius, ((@$ljsigma[@$parents[$i]])/pow(2,1./6.) -@$bondlengths[$i])*pow(2,1./6.));
    }
    
    my $index = maxIndex(@radius);
    return ($radius[$index],0.1*@$ljepsilon[@$parents[$index]]);
}

#returns the index of the maximum value in an array
sub maxIndex {
    my @a = @_;
    my $index = 0;
    for(my $i = 1; $i<=$#a; $i++){
	if($a[$i] > $a[$index]){
	    $index=$i;
	}
    }
    return $index;
}

"""



REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap',
                        'rism1d']


RUNLEAP = """source leaprc.gaff
loadamberprep {name}.prepin
check MOL
loadamberparams {name}.frcmod
SaveAmberParm MOL {name}.parm7 {name}.rst7
SavePdb MOL {name}.pdb
quit
"""


def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Runs 1D-RISM calculation
                                     .""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file. Must be in pdb format
                        acceptable by Antechamber. Must have pdb
                        extension.""")
    #Optional args
    parser.add_argument('-t', '--type',
                        help="""Solvent type""", default="OTH")
    parser.add_argument('--clean_up',
                        help=""" How should auxiliary files be treated:
                        0 - delete nothing;
                        1 - delete some [default];
                        2 - delete all but input, results, and log.
                        """, default=1, type=int)                        
    return parser.parse_args(argv)


def create_prmtop_file(name):
    """Generate topology file.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension
    Returns
    -------
    out: tuple
        Returns (logfile, prmtop_name).
        logfile is a writable file object containing std. out. prmtop_name
        is the name of prepared topology file
    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'wb')
    logfile.write(str(datetime.datetime.now()))     # timestamp
    #Firstly we use antechamber to recognize atom and bonding types, and
    #generate topology
    ante_out = subprocess.check_output(['antechamber',
                     '-i', '{}.pdb'.format(no_p_name),
                     '-fi', 'pdb',
                     '-o', '{}.prepin'.format(no_p_name), #output file
                     '-fo', 'prepi',   #output format describing each residue
                     '-c', 'bcc',      #charge method  (AM1-BCC)
                     '-s', '2',    #status info ; 2 means verbose
                     '-nc', '0',   #Net molecule charge
                     '-m', '1'],   #Multiplicity
                     cwd=p)
    logfile.write(ante_out)
    #Run parmchk to generate missing gaff force field parameters
    parm_out = subprocess.check_output(['parmchk2',
                     '-i', '{}.prepin'.format(no_p_name),
                     '-f', 'prepi',
                     '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                     cwd=p)
    logfile.write(parm_out)
    logfile.flush()
    #Run tleap to generate topology and coordinates for the molecule
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'wb') as f:
        f.write(RUNLEAP.format(name=no_p_name))
    leap_out = subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    logfile.write(leap_out)
    logfile.flush()
    prmtop_name = '{}.prmtop'.format(no_p_name)
    return logfile, prmtop_name
    
    
def generate_mdl_file(name, solvent_type):
    """Generate solvent mdl file.
    Writes modified perl script and executes it.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension
    solvent_type : string
        Either OTH, WAT or ION
    """
    custom_text = PREPDAT_PL_SCRIPT_CUSTOM_PART.format(name=name, 
                                           solvent_type=solvent_type)
    script_text = PREPDAT_PL_SCRIPT_BEGIN + custom_text + PREPDAT_PL_SCRIPT_END                                          
    with open('perldat.pl', 'wb') as f:
        f.write(script_text)
    subprocess.call(['perl', '-w', './perldat.pl'])
    

def clean_up(name, level):
    """Delete junk.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    level : {0, 1, 2}
        0 - delete nothing; 
        1 - delete ANTECHAMBER*, all water but .sh and .therm,
        .frcmod, .prepin, NEWPDB.PDB, PREP.INF, ATOMTYPE.INF, runleap.in
        sqm*, leap.log;
        2 - delete ALL but RESULTS_NAME and logfile.
    """
    p, no_p_name = os.path.split(name)
    to_del1_glob = ['ANTECHAMBER*', 'sqm*', 'water*vv*']
    to_del1_files = [no_p_name + '.prepin', no_p_name + '.frcmod',
                    'ATOMTYPE.INF', 
                    'leap.log', 'NEWPDB.PDB', 'PREP.INF', 'runleap.in']
    will_be_deleted_list = []
    if level == 1:
        for wildcard in to_del1_glob:
            will_be_deleted_list.extend(glob.glob(os.path.join(p, wildcard)))
        will_be_deleted_list.extend([os.path.join(p, f) for f in \
                                                            to_del1_files])
    if level == 2:
        all_files = os.listdir(p)
        mdl_name = '{}.mdl'.format(no_p_name)
        all_files.remove(mdl_name)
        will_be_deleted_list.extend([os.path.join(p, f) for f in all_files])
    for f in will_be_deleted_list:
        os.unlink(f)

    
def main(argv):
    args = process_command_line(argv)
    for executable in REQUIRED_EXECUTABLES:
        if not distutils.spawn.find_executable(executable):
            raise NameError("{} is not found!".format(executable))
    print('Writing and executing many scripts...')
    if args.file[-4:] == '.pdb':     
        name = args.file[:-4]
        logfile, prmtop_name = create_prmtop_file(name)
    elif args.file[-7:] == '.prmtop':
        prmtop_name = args.file
        name = prmtop_name[:-7]
    else:
        raise ValueError('Unknown file extension')
    generate_mdl_file(name, args.type)
    print('Calculation has finished')
    clean_up(name, args.clean_up)

#    clean_up(name, args.temperature, args.clean_up)



if __name__ == '__main__':
    main(sys.argv[1:])



    
    
    
    
#
#
#    #Generate water susceptibility file
#    xvv_script_name_no_p = 'water_{}_script.sh'.format(T)
#    xvv_script_name = os.path.join(p, xvv_script_name_no_p)
#    diel = round(water_dielectric_const(T), 3)
#    conc = round(water_concentration(T), 3)
#    succ_srcirpt = SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc,
#                                            wmodel=wmodel, rism1d=rism1d,
#                                            closure=closure)
#    with open(xvv_script_name, 'wb') as f:
#        f.write(succ_srcirpt)
#    xvv_out = subprocess.check_output(['bash', xvv_script_name_no_p], cwd=p)
#    logfile.write(xvv_out)
#    logfile.flush()
#    #Get compressibility from 1D-RISM therm. output
#    water_therm_name = 'water_{}.therm'.format(T)
#    water_therm_p = os.path.join(p, water_therm_name)
#    with open(water_therm_p, 'rb') as f:
#        w_lines = f.readlines()
#    compres = float(w_lines[2].split()[-1])
#    units = w_lines[2].split()[-2]
#    if units == '[1/kPa]':  # # !! In ambertools 14 compressiblity has [1/kpa] units!
#        compres = compres*1.0e4

