#!/bin/sh

DEBUG="gdb -tui --args "
DEBUG=""

$DEBUG ../bin/qelikediff -f default_nubar_qe_like_scint.txt \
    -o default_nubar_qe_like_scint.root \
    -m 20000 \
    -s -14
