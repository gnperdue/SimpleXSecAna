#!/bin/sh

DEBUG="gdb -tui --args "
DEBUG=""

nubarfilelistlist="
default_nubar_qe_like_scint.txt
zexp_nubar_qe_like_scint.txt
"

nufilelistlist="
default_nu_qe_like_scint.txt
zexp_nu_qe_like_scint.txt
default_nu_qe_like_carbon.txt
zexp_nu_qe_like_carbon.txt
default_nu_qe_like_iron.txt
zexp_nu_qe_like_iron.txt
default_nu_qe_like_lead.txt
zexp_nu_qe_like_lead.txt
"


nubarfilelistlist="
default_nubar_qe_like_scint.txt
"
nufilelistlist="
default_nu_qe_like_scint.txt
"

NEVT=-1
NEVT=20000


for filelist in $nubarfilelistlist
do
    fileroot=$(echo $filelist | cut -d. -f1)
    listname=${fileroot}.txt
    rootname=${fileroot}.root
    echo $listname
    echo $rootname
    echo ""
    $DEBUG ../bin/qelikediff -f $listname \
        -o $rootname \
        -m $NEVT \
        -s -14
done

for filelist in $nufilelistlist
do
    fileroot=$(echo $filelist | cut -d. -f1)
    listname=${fileroot}.txt
    rootname=${fileroot}.root
    echo $listname
    echo $rootname
    echo ""
    $DEBUG ../bin/qelikediff -f $listname \
        -o $rootname \
        -m $NEVT \
        -s 14
done
