#!bin/bash
for temp in $*; do
    link="$link /glusterfs/bionimbus/Human_ENCODE/$temp"
done
echo $link

echo "$(tput setaf 4 ; tput bold)F GFP seq, Left of Match$(tput sgr0)"
zgrep --color=always GATTATGATATTCCAACTA /glusterfs/bionimbus/Human_ENCODE/$link
echo "$(tput setaf 4 ; tput bold)R GFP seq, Right of Match$(tput sgr0)"
zgrep --color=always AGTTGGAATATCATAATC /glusterfs/bionimbus/Human_ENCODE/$link
echo "$(tput setaf 4 ; tput bold)F Kan seq, Right of Match$(tput sgr0)"
zgrep --color=always TTGACGAGTTCTTCTGA /glusterfs/bionimbus/Human_ENCODE/$link
echo "$(tput setaf 4 ; tput bold)R Kan seq, Left of Match$(tput sgr0)"
zgrep --color=always TCAGAAGAACTCGTCAA /glusterfs/bionimbus/Human_ENCODE/$link
