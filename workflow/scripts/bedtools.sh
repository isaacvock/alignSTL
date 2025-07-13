#!/bin/bash

input=$1
output=$2

bedtools genomecov -ibam "$input" -bg -5 -strand + | awk '{ printf "%s \t %d \t %d \t %d\n", $1,$2,$3,$4 }' > "$output"