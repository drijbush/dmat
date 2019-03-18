#!/bin/bash
# Generate DLPOLY ford-doc
# Clean old Ford Doc
rm -f *.html
rm -rf tmp

mkdir tmp
files=`ls *.f90`
for j in $files; do 
    ./ford-docgen.pl $j > tmp/${j##*/}; 
done
echo "Generating ford docs..."
ford --debug ootest.md
if [ $? -ne 0 ]; then
    echo "Error in Ford"
    exit 1
fi
echo "Done, cleaning up..."
rm -rf tmp
