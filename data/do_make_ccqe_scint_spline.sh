#!/bin/sh

# Optionally supply an extra tag for the file name.

XMLOUT=ccqe_hydrogen_carbon_splines
if [ $# -gt 0 ]; then
  XMLOUT=${XMLOUT}_$1
fi
XMLOUT=${XMLOUT}.xml
echo "Making xml file $XMLOUT"

nice gmkspl -p 14,-14 -t 1000060120,1000010010 -o ${XMLOUT} --event-generator-list CCQE


