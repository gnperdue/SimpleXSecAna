#!/bin/sh

DEBUG="gdb -tui --args "
DEBUG=""

# nubarfilelistlist="
# default_nubar_qe_like_scint.txt
# zexp_nubar_qe_like_scint.txt
# "

# nufilelistlist="
# default_nu_qe_like_scint.txt
# zexp_nu_qe_like_scint.txt
# default_nu_qe_like_carbon.txt
# zexp_nu_qe_like_carbon.txt
# default_nu_qe_like_iron.txt
# zexp_nu_qe_like_iron.txt
# default_nu_qe_like_lead.txt
# zexp_nu_qe_like_lead.txt
# "

NEVT=100000
NEVT=-1


$DEBUG ../bin/qelikediff -f default_nubar_qe_like_scint.txt \
    -o default_nubar_qe_like_scint.root \
    -m $NEVT \
    -s -14 \
    -n 0 -x 50 -c 13 -u 13
    # -n 1.5 -x 10 -c 7 -u 13

$DEBUG ../bin/qelikediff -f zexp_nubar_qe_like_scint.txt \
    -o zexp_nubar_qe_like_scint.root \
    -m $NEVT \
    -s -14 \
    -n 0 -x 50 -c 13 -u 13

$DEBUG ../bin/qelikediff -f default_nu_qe_like_scint.txt \
    -o default_nu_qe_like_scint.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 13 -u 13

$DEBUG ../bin/qelikediff -f default_nu_qe_like_carbon.txt \
    -o default_nu_qe_like_carbon.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 12 -u 12

$DEBUG ../bin/qelikediff -f default_nu_qe_like_iron.txt \
    -o default_nu_qe_like_iron.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 56 -u 56

$DEBUG ../bin/qelikediff -f default_nu_qe_like_lead.txt \
    -o default_nu_qe_like_lead.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 207 -u 207

$DEBUG ../bin/qelikediff -f zexp_nu_qe_like_scint.txt \
    -o zexp_nu_qe_like_scint.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 13 -u 13

$DEBUG ../bin/qelikediff -f zexp_nu_qe_like_carbon.txt \
    -o zexp_nu_qe_like_carbon.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 12 -u 12

$DEBUG ../bin/qelikediff -f zexp_nu_qe_like_iron.txt \
    -o zexp_nu_qe_like_iron.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 56 -u 56

$DEBUG ../bin/qelikediff -f zexp_nu_qe_like_lead.txt \
    -o zexp_nu_qe_like_lead.root \
    -m $NEVT \
    -s 14 \
    -n 0 -x 50 -c 207 -u 207

