#! /bin/sh

FILE=$1
export FILE

DX=`awk '{print NF}' $FILE | head -n1`
export DX
DY=`wc -l $FILE`
export DY

idl << EOF

.r plot_mask.pro
plot_mask, 1
EOF

gv field.ps &
