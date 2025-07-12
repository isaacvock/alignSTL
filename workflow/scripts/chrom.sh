input=$1
output=$2

samtools view -H "$input" \
    | awk -v OFS="\t" ' $1 ~ /^@SQ/ {split($2, chr, ":")
                                        split($3, size, ":")
                                        print chr[2], size[2]}' > "$output"