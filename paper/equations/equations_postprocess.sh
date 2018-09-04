#!/bin/bash -e
# cd to this directory
for f in `find . -name "*.txt"`; do
  # Zeilenumbruch ersetzen
  sed -i ':a;N;$!ba;s/\n/ /g' $f
  # Latex-Ausdrücke für Darstellung optimieren
  # it unnätig
  sed -i "s/\\\\it//g" $f
  # nx, ... mit Index
  sed -i "s/ nx/n_x/g" $f
  sed -i "s/ ny/n_y/g" $f  
  sed -i "s/ nz/n_z/g" $f
  sed -i "s/ ox/o_x/g" $f
  sed -i "s/ oy/o_y/g" $f  
  sed -i "s/ oz/o_z/g" $f
  sed -i "s/ ax/a_x/g" $f
  sed -i "s/ ay/a_y/g" $f  
  sed -i "s/ az/a_z/g" $f
  # Cos/Sin abkürzen mit Index
  sed -i "s/S/S_/g" $f
  sed -i "s/C/C_/g" $f
  sed -i "s/delta/\{\\\\delta\}/g" $f
done
