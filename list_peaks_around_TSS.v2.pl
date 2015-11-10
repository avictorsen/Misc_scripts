#!/usr/bin/perl
use strict;

#Makes hash of good data files containing region peaks

my %hash;
my @gooddata;
#number of bases on either side:
my $n = 2500;
my $step = 100;

my $homedir = "/glusterfs/data/modencode/WhiteLab_UniformProcessing/hg/processing/";
Search($homedir);

#make array of all good bed files
sub Search{
	opendir (DIR, $homedir) || die;
	my @files = grep { !/^\.{1,2}/ } readdir(DIR);
	closedir (DIR);
#       foreach (@files){print "$_\n";}
        for (@files){
                if (-d $homedir."/".$_){
                        #print "$_\n";
                        my $batch = $_;
                        my $newdir = $homedir."/".$batch."/";
                        my $newdir= $newdir."wig/";#change if you want to search other sub$
                        opendir (SUBDIR, $newdir);#||die $newdir;

                        #get all files that match pattern
                        my @files2 = grep {/^\S+Rep0\.tagAlign_VS\S+Input_Rep0\.tagAlign.w$
                        #my @files2 = grep {/^\S+Rep0\.tagAlign_density.wig/} readdir(SUBD$
                        #my @files2 = grep {/spp\.optimal.*IP_Rep0.*Input.*bed$/} readdir $
                        #foreach (@files2){print "$_\n";}<STDIN>;

                        foreach (@files2){
                                my $temp = $_;
                                #print $temp;<STDIN>;
                                my $gene;
                                #if ($temp =~ /^spp\.optimal\.(.*?)-/){$gene = $1;}
                                if ($temp =~ /^(.*?)_/){$gene = $1;}
                                $temp = $newdir.$temp;
                                foreach my $item (@gooddata){
                                        chomp $item;
                                        if (lc $item eq lc "$gene $batch"){
                                                #print "X $item X\tY $gene $batch Y\t";
                                                for (my $i = 0; $i < scalar @gooddata; $i+$
                                                        if (lc $gooddata[$i] eq lc $item){
                                                                my $temp = splice (@goodda$
                                                        }
                                                }
                                                #print "Z $temp Z\n";
                                                $hash{$item} = $temp;
                                        }
                                #else{ print "X $item X\tY $gene $batch Y\n";}
                                }
                                undef $temp;undef $gene;
                        }
                        closedir SUBDIR;undef $batch; undef $newdir;
                }
        }
}

print "could not find files for:\n";
foreach (@gooddata){print "$_\n";}
print "hit return";<STDIN>;
#print "\nmatched files:\n";
#foreach (keys %hash){print "$_\t$hash{$_}\n";}
#<STDIN>;



#foreach (keys %hash){print "$_\t$hash{$_}\n";}
#die;

open (OUTPUT, ">list_peaks_around_TSS_signal_output.txt")||die;
open (OUTPUT2, ">list_peaks_around_TSS_background_subtracted_signal_output.txt")||die;
open (OUTPUT3, ">list_peaks_around_TSS_background_output.txt")||die;
open (OUTPUT4, ">list_peaks_around_TSS_coefficient_of_variation_output.txt")||die;

print OUTPUT "\t\t";
print OUTPUT2 "\t\t";
print OUTPUT3 "\t\t";

#print headers in output files
for (my $i = -$n; $i < $n; $i += $step){print OUTPUT "$i\t"; print OUTPUT2 "$i\t"; print OUTPUT3 "$i\t";}
 
#get genome lengths
#these genome lenghts are based on :/raid/database/ftp.flybase.net/releases/FB2011_09/dmel_r5.41/fasta/dmel5.41/dmel-all-chromosome-r5.41.fasta

my %chromosomes = (
		'chr1'=> '0',
		'chr2'=> '249250621',
		'chr3'=> '492449994',
		'chr4'=> '690472424',
		'chr5'=> '881626700',
		'chr6'=> '1062541960',
		'chr7'=> '1233657027',
		'chr8'=> '1392795690',
		'chr9'=> '1539159712',
		'chr10'=> '1680373143',
		'chr11'=> '1815907890',
		'chr12'=> '1950914406',
		'chr13'=> '2084766301',
		'chr14'=> '2199936179',
		'chr15'=> '2307285719',
		'chr16'=> '2409817111',
		'chr17'=> '2500171864',
		'chr18'=> '2581367074',
		'chr20'=> '2659444322',
		'chr19'=> '2722469842',
		'chr22'=> '2781598825',
		'chr21'=> '2832903391',
		'chrX'=> '2881033286',
		'chrY'=> '3036303846');
    
my @geneTSSsarray; 
#open file with TSSs
open (GENES, "<./refGene.txt")||die;

while (<GENES>){
	#pull attrbutes for each gene
        my $chr; my $strand; my $start; my $end;
	my @line = split("\t", $_);
	if ($line[2] =~ /^chr/){$chr = $line[2];}
	if ($line[3] =~ /(.+?)/){$strand = $1;}
	if ($line[4] =~ /(\d+)/){$start = $1;}
	if ($line[5] =~ /(\d+)/){$end = $1;}
	if ($strand eq "-"){$start = $end;}
	my $point = int(($chromosomes{$chr} + $start)/$step)*$step;
	#print $point;<STDIN>;
	if ($point != 0){push (@geneTSSsarray, $point);}
	#print "$chr \t $start \t $end \t $strand \t $point \n------\n\n";<STDIN>;
}
close GENES;
die;

my $m = 0;
#open wig file for each good data set
foreach (keys %hash){
	$m++;
	#get gene-ish name from file
	my $filename = $hash{$_};
	$filename  =~ /\/([\d\w\-]+?)-[GM]/;
	my $genename = $1;
	open (READ, "<".$hash{$_})||die;
	my %wighash;
	#make hash of coordinates
	my @wigfile = <READ>;
	shift @wigfile;
	print "\n\nprocessing wig $m of ".(keys %hash)." :\n";
	foreach (@wigfile){
                my @line = split(" ", $_);
		if ($chromosomes{$line[0]}){
			#print "$line[0] \t $line[1]\r";
			my $newcoor = int(($chromosomes{$line[0]} + $line[1])/100) * 100;
			$wighash{$newcoor} = $line[3];
			#print "$line[0] \t $line[1]\t :: \t $chromosomes{$line[0]}\t$newcoor\n";<STDIN>; 
		}
		#else{ print "$line[0] \t $line[1]\r";}
	}

        #calculate random wig enrichment
        #my %randohash;
        #print "creating random hash\n";
        #for (my $i = -$n; $i < $n; $i += $step){
        #        $randohash{$i} = 0;
        #}
        #my $j = 0;
	#my @randosample;
        #for (my $i = 0; $i < scalar @geneTSSsarray; $i++){
        #        #print "$j of ".scalar @geneTSSsarray."\r";
        #        my $base = int(rand($chromosomes{'end'})/100)*100;
        #        for (my $i = -$n; $i < $n; $i += $step){
        #                if ($i == 0){push (@randosample, $wighash{$base});}
	#		$randohash{$i} = $randohash{$i} + $wighash{$base + $i};
        #        }
        #        $j++;
        #}

	#calculating coefficient of variation for dataset
	#my $CVsum = 0;
	#foreach (@randosample){$CVsum += $_;}
	#my $CVmean = $CVsum / scalar @randosample;
	#my $CVstd;
	#foreach (@randosample){$CVstd += ( $_ - $CVmean ) **2;}
	#$CVstd = ( $CVstd / scalar @randosample ) ** 0.5;
	#my $CV = $CVstd / $CVmean;
	#print OUTPUT4 "$genename\t$filename\t$CVmean\t$CVstd\t$CV\n";
	#my $bigavg = 0;
	#my $littleavg = 0;
	#foreach (keys %randohash){
	#	$littleavg += $randohash{$_}/scalar @geneTSSsarray;
	#	$bigavg += $randohash{$_}
	#}
	#$bigavg = $bigavg / (( $n / $step ) * 2);
	#$littleavg = $littleavg / (( $n / $step ) * 2);
	#my $std = 0;
	#foreach (keys %randohash){
	#	$std += (($randohash{$_}/scalar @geneTSSsarray - $littleavg)**2);
	#}
	#$std = ( $std / (( $n / $step ) * 2))**(0.5);
	#my $min = $littleavg + 2 * $std;


	##reset TSS hash to 0's
	##print "\nresetting hash\n";
	my $tempn = $step;
	my $background;
	my %TSShash;
        for (my $i = -$n; $i < $n; $i += $step){
  	      $TSShash{$i} = 0;
        }

 	do{
		$tempn += $step;
		print "checking bases $tempn around TSSs for $genename.";
		#check every TSS
		my $j = 0;
		$background = 0;
		foreach (@geneTSSsarray){
			#print "$j of ".scalar @geneTSSsarray."\r";
			my $point = int($_/100)*100;
			$background += $wighash{$point - 10000}; 
			#print "\n$point\t$_";<STDIN>;
			for (my $i = -$tempn; $i < $tempn; $i += $step){
				#print $point+$i;<STDIN>;
				if ($wighash{$point +$i}){$TSShash{$i} = $TSShash{$i} + $wighash{$point + $i};}
			}
			$j++;
		print "\r";
		}
	}while($tempn < 2500); # && $TSShash{-$tempn} > $background);
	print OUTPUT "\n$genename\t$filename\t";
        for (my $i = -$n; $i < $n; $i += $step){print OUTPUT ($TSShash{$i}."\t");}

        #print OUTPUT2 "\n$genename\t$filename\t";
        #for (my $i = -$n; $i < $n; $i += $step){print OUTPUT2 ($TSShash{$i}-$bigavg)."\t";}

	#print OUTPUT3 "\n$genename\t$filename\t";
	#for (my $i = -$n; $i < $n; $i += $step){print OUTPUT3 ($randohash{$i}."\t");}

}
