#!/bin/sh

if [ -z $XSECSPLINEDIR ]; then
  echo "No XSECSPLINEDIR environment variable defined! Cannot proceed!"
  exit 0
fi

if [ ! -d $XSECSPLINEDIR ]; then
  echo "No XSECSPLINEDIR directory exists! Cannot proceed!"
  exit 0
fi

cp ccqe_hydrogen_carbon_splines.xml $XSECSPLINEDIR
