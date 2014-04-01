#!/bin/sh

# Generic run on oxygen to test the installation.

NUMEVT=100
if [ $# -gt 0 ]; then
  NUMEVT=$1
fi

if [ ! -f $XSECSPLINEDIR/gxspl-vA-v2.8.0.xml ]; then
  echo "Cross section spline file is missing."
  exit 0
fi

gevgen -n $NUMEVT -p 14 -t 1000080160 -e 0,10 -r 100 \
  -f 'x*exp(-x)' \
  --seed 2989819 --cross-sections $XSECSPLINEDIR/gxspl-vA-v2.8.0.xml >& run_log.txt

