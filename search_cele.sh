#!bin/bash
for temp in $*; do
    link="$link /glusterfs/bionimbus/modENCODE_ChIP-seq/$temp"
done
echo $link

echo "$(tput setaf 4 ; tput bold)F GFP seq, Left of Match$(tput sgr0)"
zgrep --color=always GAAGTGCATACCAATCA /glusterfs/bionimbus/modENCODE_ChIP-seq/$link
echo "$(tput setaf 4 ; tput bold)R GFP seq, Right of Match$(tput sgr0)"
zgrep --color=always TGATTGGTATGCACTTC /glusterfs/bionimbus/modENCODE_ChIP-seq/$link
echo "$(tput setaf 4 ; tput bold)F Kan seq, Right of Match$(tput sgr0)"
zgrep --color=always AAGGATGACGACGACAAG /glusterfs/bionimbus/modENCODE_ChIP-seq/$link
echo "$(tput setaf 4 ; tput bold)R Kan seq, Left of Match$(tput sgr0)"
zgrep --color=always CTTGTCGTCGTCATCCTT /glusterfs/bionimbus/modENCODE_ChIP-seq/$link
