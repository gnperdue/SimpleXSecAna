#!/bin/sh

DEBUG=""
DEBUG="gdb -tui --args "

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

NEVT=-1
NEVT=100000


$DEBUG ../bin/qelikediff -f default_nubar_qe_like_scint.txt \
    -o default_nubar_qe_like_scint.root \
    -m $NEVT \
    -s -14 \
    -n 1.5 -x 10 -c 7
    # -n 0 -x 50 -c 13 

# $DEBUG ../bin/qelikediff -f default_nu_qe_like_scint.txt \
#     -o default_nu_qe_like_scint.root \
#     -m $NEVT \
#     -s 14 \
#     -n 1.5 -x 10 -c 6
