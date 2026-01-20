#!/bin/bash

# Downloads metadata file and replaces white spaces w/ underscores

link="$1"
out_name="$2"

#*************************
# DOWNLOAD METADATA FILE *
#*************************

wget "$link" -O "$out_name"
sed -e 's/ /\_/g' "$out_name" > tmp
mv tmp "$out_name"
