#!/bin/bash
if [ ! -f "OUTCAR" ];
then
echo "OUTCAR file not present"
fi
if [ -f "OUTCAR" ];
then
echo "OUTCAR file present in current directory"
fi
if [ ! -f "pos-force.txt" ];
then
echo "pos-force.txt file not present, execute reduce.x to get this file"
fi
if [ -f "pos-force-new.txt" ];
then
rm pos-force-new.txt
fi
if [ -f "pos-force.txt" ];
then
echo "pos-force.txt file present in the current directory"
echo " # POSITION      TOTAL FORCE" >> pos-force-new.txt
echo "1       `cat OUTCAR |grep "free  energy   TOTEN  ="|awk -F "=" '{print $2}'` =t, E(t)" >> pos-force-new.txt
cat pos-force.txt >> pos-force-new.txt
mv pos-force-new.txt pos-force.txt
fi 
