#!/usr/bin/perl -w 
use strict ;
use warnings ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::File ;

opendir (DIR, "/glusterfs/bionimbus/modENCODE_ChIP-seq/");
open (OUTPUT, ">output.txt");
print OUTPUT "BID\tMACHINE\tFLOWCELL\tLANE\tINDEX\n";

my @IDs = @ARGV;
my @files = readdir(DIR);

foreach (@IDs){
	my $BID = $_;
	my $query = $BID."_";
	foreach (@files){
		if ($_ =~ m/^$query/){
                        my @tempname = split('_', $_);
			my $z = new IO::Uncompress::Gunzip "/glusterfs/bionimbus/modENCODE_ChIP-seq/".$_ or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
			my $line = $z->getline();
			#print $line;
			my @array = split (':', $line);
			if ($array[9] eq '0'){
				print OUTPUT "$tempname[2]\t$BID\t$array[0]\t$array[2]\t$array[3]\t$array[10]";
			}
			else{
				print OUTPUT "$tempname[2]\t$BID\t$array[0]\t$array[2]\t$array[3]\t$array[9]";
			}
		}
	}
}
