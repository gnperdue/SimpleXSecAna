#!/bin/sh

# Generic run on oxygen to test the installation.

NUMEVT=100000
if [ $# -gt 0 ]; then
  NUMEVT=$1
fi

if [ ! -f $XSECSPLINEDIR/gxspl-vA-v2.8.0.xml ]; then
  echo "Cross section spline file is missing."
  exit 0
fi

gevgen -n $NUMEVT -p 14 -t 1000060120 -e 3 -r 51 \
  --seed 2989819 --cross-sections $XSECSPLINEDIR/gxspl-vA-v2.8.0.xml \
  --event-generator-list CCQE >& run_log.txt

