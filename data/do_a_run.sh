#!/bin/sh

if [ -z $XSECSPLINEDIR ]; then
  echo "No XSECSPLINEDIR environment variable defined! Cannot proceed!"
  exit 0
fi

if [ ! -d $XSECSPLINEDIR ]; then
  echo "No XSECSPLINEDIR directory exists! Cannot proceed!"
  exit 0
fi

NUMEVT=1000
if [ $# -gt 0 ]; then
  NUMEVT=$1
fi

SPLINE="$XSECSPLINEDIR/ccqe_hydrogen_carbon_splines.xml"

if [ ! -f $SPLINE ]; then
  echo "Cross section spline file is missing."
  exit 0
fi

gevgen -n $NUMEVT -p 14 -t 1000060120,1000010010 -e 0,10 -r 100 \
  -f 'x*exp(-x)' \
  --seed 2989819 \
  --cross-sections $SPLINE \
  --event-generator-list CCQE >& run_log.txt

