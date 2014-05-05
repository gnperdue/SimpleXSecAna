#!/bin/sh

FIRST=100
LAST=100
for i in $(eval echo "{$FIRST..$LAST}")
do
  echo "$i"
  gntpc -i ./gntp.$i.ghep.root -f rootracker
  gntpc -i ./gntp.$i.ghep.root -f gst
done

ls -1 *gtrac.root > gtrac_list.txt
ls -1 *gst.root > gst_list.txt

if test -d ../list; then echo ListOk >& /dev/null ; else mkdir -p ../list; fi

mv gtrac_list.txt ../list
mv gst_list.txt ../list
