#!/bin/bash

# split the master vcf file into subsets of 100k each
split -l 100000 /genomes/species/spret-b6/SPRET_EiJ.variants/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf ./SPRET_EiJ.snp.100k.subset_

# get vcf header
head - n 69 ./SPRET_EiJ.snp.100k.subset_aa > snp.header.txt

# Add header to each file and save as vcf
# Loop through the range of aa to rl
for letter in {a..r}; do
    for num in {a..l}; do
        filename="SPRET_EiJ.snp.100k.subset_${letter}${num}"
        output_file="${filename}.vcf"

        # Check if the file exists before appending the header
        if [ -f "$filename" ]; then
            cat snp.header.txt "$filename" > "$output_file"
            echo "Appended header to $filename and saved as $output_file"
        else
            echo "File $filename does not exist."
        fi
    done
done

