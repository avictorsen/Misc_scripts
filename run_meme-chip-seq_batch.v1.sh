#!bin/bash
bindir="/mnt/cinder/SCRATCH"
mkdir HOT
mkdir COLD
mkdir fasta
mkdir controls
mkdir masked

for file in $1/*.bed; do
    filename="${file##*/}"
    echo "$(tput setaf 4 ; tput bold)Processing $filename$(tput sgr0)"

    #fix bedfiles
    echo "$(tput setaf 4 ; tput bold)Fixing bedfiles$(tput sgr0)"
    cat $file | awk '{
        FS="\t";
        if ($3 ~ /e/){
            $3 = sprintf("%d",$3);
        };
        if ($2 ~ /e/){
            $2 = sprintf("%d",$2);
        };
        if ($1 ~ /dmel_mito/){
            $1 = "M";
        }
        if ($1 ~ /Uextra/){
            $1 = "UEXTRA";
        }
        if ($1 ~ /Het/){
            gsub(/Het/,"HET",$1);
        }
        else{;}
        start=$2;
        end=$3;
        if (start < end){
            temp=sprintf("%.0f", (end - start) / 2);
            $3=temp + 250 + start;
            $2=temp - 250 + start;
        }
        else{
            temp=sprintf("%.0f", (start - end) / 2);
            $2=temp - 250 + end;
            $3=temp + 250 + end;
        }
        print $0;
        }' OFS="\t" > ./temp.bed

    bedtools sort -i ./temp.bed > ./temp2.bed
    bedtools merge -d 0 -i ./temp2.bed > ./temp3.bed

    #make cold files
    #echo "$(tput setaf 4 ; tput bold)Removing HOTspots$(tput sgr0)"
    bedtools slop -i ./temp3.bed -g $bindir/dm3.genome -b 0 >> ./HOT/$filename
    #bedtools intersect -v -a ./HOT/$filename -b $bindir/HOTregions.bed > ./COLD/$filename.cold
    rm -f ./temp*.bed

    #make fasta files
    echo "$(tput setaf 4 ; tput bold)Writing fasta files$(tput sgr0)"
    #bedtools getfasta -fi $bindir/dm3.fa -bed ./COLD/$filename.cold -fo ./fasta/$filename.cold.fa
    bedtools getfasta -fi $bindir/dm3.fa -bed ./HOT/$filename -fo ./fasta/$filename.fa

    #make controls
    echo "$(tput setaf 4 ; tput bold)Generating controls$(tput sgr0)"
    #bedtools shuffle -chrom -incl $bindir/Peaks.small.bed -i ./COLD/$filename.cold -g $bindir/dm3.genome > ./temp.bed
    #bedtools sort -i ./temp.bed > temp2.bed
    #bedtools merge -d 0 -i temp2.bed > ./controls/$filename.cold.shuf
    #rm -f ./temp*.bed
    #bedtools getfasta -fi $bindir/dm3.fa -bed ./controls/$filename.cold.shuf -fo ./controls/$filename.cold.shuf.fa

    bedtools shuffle -chrom -incl $bindir/Peaks.small.bed -i ./HOT/$filename -g $bindir/dm3.genome > ./temp.bed
    bedtools sort -i ./temp.bed > temp2.bed
    bedtools merge -d 0 -i temp2.bed > ./controls/$filename.shuf
    rm -f ./temp*.bed
    bedtools getfasta -fi $bindir/dm3.fa -bed ./controls/$filename.shuf -fo ./controls/$filename.shuf.fa

    #mask all files
    echo "$(tput setaf 4 ; tput bold)Masking tandem repeats$(tput sgr0)"
    echo "$(tput setaf 1 ; tput bold)Masking all-peaks file$(tput sgr0)"
    ~/TOOLS/RepeatMasker/RepeatMasker ./fasta/$filename.fa -species drosophila -pa 8 -dir ./masked -s
    echo "$(tput setaf 1 ; tput bold)Masking all-peaks control file$(tput sgr0)"
    ~/TOOLS/RepeatMasker/RepeatMasker ./controls/$filename.shuf.fa -species drosophila -pa 8 -dir ./masked -s
    #echo "$(tput setaf 1 ; tput bold)Masking cold peaks file$(tput sgr0)"
    #~/TOOLS/RepeatMasker/RepeatMasker ./fasta/$filename.cold.fa -species drosophila -pa 8 -dir ./masked -s
    #echo "$(tput setaf 1 ; tput bold)Masking cold peaks control file$(tput sgr0)"
    #~/TOOLS/RepeatMasker/RepeatMasker ./controls/$filename.cold.shuf.fa -species drosophila -pa 8 -dir ./masked -s

    #run ame
    #echo "$(tput setaf 4 ; tput bold)Running ame$(tput sgr0)"
    #TF=`echo $filename | sed 's/spp.optimal.//' | grep -P -o "^\S+?_WA"`
    #TF=${TF%_WA}
    #echo "$(tput setaf 1 ; tput bold)processing all-peaks$(tput sgr0) $TF"
    #ame --o ./$TF --control ./masked/$filename.shuf.fa.masked ./masked/$filename.fa.masked ~/TOOLS/meme/motif_databases/JASPAR/JASPAR_CORE_2014_insects.meme
    #TF=$TF"_cold"
    #echo "$(tput setaf 1 ; tput bold)processing cold-peaks$(tput sgr0) $TF"
    #ame --o ./$TF --control ./masked/$filename.cold.shuf.fa.masked ./masked/$filename.cold.fa.masked ~/TOOLS/meme/motif_databases/JASPAR/JASPAR_CORE_2014_insects.meme

    #run meme
    #echo "$(tput setaf 4 ; tput bold)Running meme on all peaks$(tput sgr0):"
    #TF=`echo $filename | sed 's/spp.optimal.//' | grep -P -o "^\S+?_WA"`
    #TF=${TF%_WA}
    #meme ./masked/$filename.fa.masked -dna -p 8 -oc ./$TF -mod zoops -evt 0.05 -maxsize 10000000

    #run centrimo
    #CTF=$TF"_meme-Cent"
    #centrimo --local --oc ./$CTF --neg ./masked/$filename.shuf.fa.masked ./masked/$filename.fa.masked ./$TF/meme.txt


    #run meme-chip
    echo "$(tput setaf 4 ; tput bold)Running MEME-ChIP on all peaks$(tput sgr0):"
    TF=`echo $filename | sed 's/spp.optimal.//' | grep -P -o "^\S+?_WA"`
    TF=${TF%_WA}
    ATF=$TF"_meme-chip_all"
    nmeme=`grep ">" ./masked/$filename.shuf.fa.masked | wc -l`
    meme-chip -oc ./$ATF -meme-p 8 -meme-maxsize 10000000 -meme-mod zoops -group-thresh 0.001 -filter-thresh 0.0000000001 -bfile $bindir/HOTregions.HMM.background -neg ./masked/$filename.shuf.fa.masked ./masked/$filename.fa.masked

    #echo "$(tput setaf 4 ; tput bold)Running MEME-ChIP on cold peaks$(tput sgr0):"
    #CTF=$TF"_meme-chip_cold"
    #nmeme=`grep ">" ./masked/$filename.cold.shuf.fa.masked | wc -l`
    #meme-chip -oc ./$CTF -meme-p 8 -meme-maxsize 10000000 -nmeme $nmeme -meme-mod zoops -norand -neg ./masked/$filename.cold.shuf.fa.masked ./masked/$filename.cold.fa.masked

    echo "$(tput setaf 4 ; tput bold)--------------------------------------------------------------------------------\n$(tput sgr0)"
done
