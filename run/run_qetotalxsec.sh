#!/bin/sh

DEBUG="gdb -tui --args "
DEBUG=""

NEVT=-1
NEVT=50000

$DEBUG ../bin/qetotalxsec -f default_nubar_qe_like_scint.txt \
    -o totalxsec_default_nubar_qe_like_scint.root \
    -m $NEVT

$DEBUG ../bin/qetotalxsec -f default_nu_qe_like_scint.txt \
    -o totalxsec_default_nu_qe_like_scint.root \
    -m $NEVT

$DEBUG ../bin/qetotalxsec -f default_nu_qe_like_carbon.txt \
    -o totalxsec_default_nu_qe_like_carbon.root \
    -m $NEVT

$DEBUG ../bin/qetotalxsec -f default_nu_qe_like_iron.txt \
    -o totalxsec_default_nu_qe_like_iron.root \
    -m $NEVT
