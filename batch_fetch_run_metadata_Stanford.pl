#!/usr/bin/perl -w 
use strict ;
use warnings ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::File ;

opendir (DIR, "/glusterfs/data/modencode/WhiteLab_UniformProcessing/ce/raw/Stanford_data");
open (OUTPUT, ">output.txt");
print OUTPUT "file\tMACHINE\tFLOWCELL\tLANE\tINDEX\n";

my @files = readdir(DIR);

foreach (@files){
	my $tempname = $_;
	my $z = new IO::Uncompress::Gunzip "/glusterfs/data/modencode/WhiteLab_UniformProcessing/ce/raw/Stanford_data/".$_ or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	my $line = $z->getline();
	#print $line;
	my @array = split (':', $line);
	if ($array[8] eq '0'){
		print OUTPUT "$tempname\t$array[0]\t$array[2]\t$array[3]\tERR\n";
	}
	else{
		print OUTPUT "$tempname\t$array[0]\t$array[2]\t$array[3]\t$array[8]\n";
	}
}


