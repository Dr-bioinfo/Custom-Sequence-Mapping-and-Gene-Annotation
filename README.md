# Custom Sequence Mapping and Gene Annotation

## Overview
This Python script is designed to facilitate custom sequence mapping and gene annotation tasks for genomic analysis. It automates the extraction and mapping of non-targeted sequences from a reference FASTA file using various bioinformatics tools, including samtools, minimap2, and bedtools. The script also performs gene annotation and generates summary reports. By providing the required input file paths and executing the script, users can seamlessly extract, map, and annotate non-targeted sequences, making it a valuable tool for genomics research and analysis.

## Research Question
This toolkit was developed to address questions like:
- How can we identify non-targeted sequences within a set of genomic reads?
- What are the affiliations of these non-targeted sequences to known genomes or regions?
- Are these non-targeted sequences aligning with specific genes or regions of interest within the reference genome?
- What is the depth of coverage for these aligned sequences in the reference genome?

## Prerequisites
- Conda environment named "mapping" with samtools, minimap2, and bedtools installed. The .yml file has been provided 

    ```bash
    $ conda activate mapping
    ```

## Command
To run the script, execute the following command:

```bash
$ python non_target_auto.py -f1 /path/arg/file -f2 /path/emu/csv -f3 /path/to/itscsv/ -sam /path/to/PC.fasta -seq /path/to/the/refseq/ -gb /path/to/the/genank/ -targ /path/to/target.txt -non /path/to/non_target.fasta -fin /path/to the /finaal_csv
```

## Customizable for Any Organism
This toolkit is customizable for use with any organism. You can specify the reference sequence and GenBank file as arguments to tailor the analysis to your specific research needs.

## Dependencies
- argparse
- sys
- csv
- subprocess
- pandas
- Bio (Biopython)

## Workflow
1. **Read ARG and 16S Reads**: The script reads ARG headers from the provided `resfinder.tsv` file and 16S headers from the `emu_best_hits.csv` & `its_best_hits.csv` file.

2. **Combine Headers**: It combines the extracted headers and writes them to a file called `targeted.txt`.

3. **Compare with Reference**: The script compares the combined headers with the reference sequences provided in the `PC.fasta` file and extracts non-matching headers and sequences. Non-targeted sequences are saved in the `non_targeted_read.fasta` file.

4. **Index Reference Genome**: The reference genome is indexed using minimap2, and the index is saved as `ref.mmi`.

5. **Map Non-Targeted Sequences**: Non-targeted sequences are mapped to the reference genome using minimap2, and the output is saved as `non_targeted.sam`.

6. **Convert SAM to BAM**: SAM files are converted to BAM format using samtools.

7. **Sort BAM**: BAM files are sorted using samtools, and the sorted output is saved as `non_targeted_sorted.bam`.

8. **Index Sorted BAM**: The sorted BAM file is indexed using samtools.

9. **Create Pileup File**: A pileup file is created from the sorted BAM file and saved as `non_targeted_count.tsv`.

10. **Convert BAM to BED**: BAM files are converted to BED format using bedtools, and the output is saved as `non_targeted_sorted.bed`.

11. **Extract Mapped Positions**: Mapped positions are extracted from the BED file and saved as `mapped_positions.csv`.

12. **Extract Gene Information**: Gene information is extracted from the provided GenBank file and mapped to the pileup file. The result is saved as `mapped_genes_info.csv`.

13. **Merge Mapped Positions and Gene Information**: Mapped positions and gene information are merged into a final output CSV file specified by the `-fin` argument.

## Output
The script generates the following files:
- `targeted.txt`: Combined header file.
- `non_targeted_read.fasta`: Non-targeted sequences.
- `ref.mmi`: Reference genome index.
- `non_targeted.sam`: SAM file.
- `non_targeted.bam`: BAM file.
- `non_targeted_sorted.bam`: Sorted BAM file.
- `non_targeted_count.tsv`: Pileup file.
- `non_targeted_sorted.bed`: BED file.
- `mapped_positions.csv`: Extracted mapped positions.
- `mapped_genes_info.csv`: Extracted gene information.
- `mapped_gene_final_output.csv`: Merged mapped positions and gene information.

##### Have a GOOD DAY 
