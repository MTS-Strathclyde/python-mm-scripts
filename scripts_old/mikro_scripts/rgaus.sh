#!/usr/bin/perl

#See on modifitseeritud Rg98 skripti versioon

# See on skriptike gaussiani toode sabasse saatmiseks
# Aktsepteerib EHK KUNAGI sisendina ka loendit. Tehtud??
chomp( $cwd = `pwd` );
if ($cwd =~ /home2/) { die "Sipsikul tööde käivitamisel peavad failid olema /home/TEIENIMI kataloogis (või selle alamkataloogis).\n"; }
while ($In = pop(@ARGV)) {
	$i = $In . ".com";
	if (-e $In || -e $i) {
		($Name) = split ".com",$In;
# Korrigeerime %NProc väärtuse!
		$i = $Name . ".com";
		$j = $Name . ".1com";
		open (IN,"$i") || die "Can't open $i: $!\n";
		open (OUT,">$j") || die "Can't open $j: $!\n";
		while (<IN>) {
			if (/nproc/i) {
				next;
			} elsif (/\#/) {
				print OUT "\%NProc=4\n";
			}
			print OUT $_;
		}
		close IN;
		close OUT;
		system "mv -f $j $i";
# Korrigeeritud!
		open (OUT,">$Name") || die "Can't open $Name: $!\n";
		print OUT "\#!/bin/bash \n\#PBS -l nodes=1:ppn=4\n";
		print OUT "\#PBS -q default\n";
	 	#print OUT "\#PBS -l walltime=10:00:00\n";
		#print OUT "\#PBS -l vmem=4gb\n";
		print OUT "cd $cwd\n/usr/local/bin/rungauss.1 $Name\n";
		print OUT "rm -f $Name\n";
		close OUT;
		system "qsub $Name";
	} else {
		print "Faili $In pole!\n"
	}
}

