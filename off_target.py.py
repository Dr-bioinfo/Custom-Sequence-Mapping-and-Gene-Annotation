## Run this code in "Conda ID_project" as it need samtools, minimap2, bed tools 
## Command : 
## $ conda activate mapping 

## $python non_target_auto.py  -targ /path/to/target_header.csv -sam /path/to/PC.fasta -seq /path/to/the/refseq/ -gb /path/to/the/genank/ 
#                              -non /path/to/output/non_target.fasta  -fin /path/to the /finaal_csv`

## Importing the Libraries 
import argparse
import sys
import csv
import subprocess
import pandas as pd
from Bio import SeqIO
import os

## Reading the targeted headers 
def read_csv_file(file_path):
    headers = []
    with open(file_path, "r") as file:
        next(file)
        for line in file:
            header = line.strip().split(",")[0]  # Extracting the 16S header from the first column
            headers.append(header)
    return headers

## Compare them with the ref fasta and extract the unique reads 
def compare_with_ref(ref_fasta_file, targeted_headers):
    ref_headers = []
    ref_sequences = []
    current_header = None
    current_sequence = ""
    with open(ref_fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    ref_headers.append(current_header)
                    ref_sequences.append(current_sequence)
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_header and current_sequence:  # For the last sequence in the file
            ref_headers.append(current_header)
            ref_sequences.append(current_sequence)

    # Remove duplicates from ref_headers and ref_sequences
    ref_headers, ref_sequences = zip(*set(zip(ref_headers, ref_sequences)))
    return ref_headers, ref_sequences

# Extract Mapped regions from BED file 
def extract_mapped_regions(bed_file):
    mapped_positions = []
    with open(bed_file, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            read_name = parts[3]
            start = int(parts[1])
            end = int(parts[2])
            mapped_positions.append({
                "Read Name": read_name,
                "Start": start,
                "End": end,
            })

    # Write mapped positions to a CSV file
    csv_output_file = "mapped_positions.csv"
    with open(csv_output_file, "w", newline='') as csvfile:
        fieldnames = ["Read Name", "Start", "End"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(mapped_positions)

    print("Mapped positions extracted and saved!")

# Extract Mapped gene 
def extract_gene_from_mapped(ref_genbank_file, pileup_out):
    mapped_genes = {}
    # Parse the genebank file
    with open(ref_genbank_file, "r") as genebank_handle:
        for record in SeqIO.parse(genebank_handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    gene_location = feature.location
                    gene_sequence = gene_location.extract(record.seq)
                    gene_name = feature.qualifiers.get("gene", ["Unknown"]) [0]
                    mapped_genes[(gene_location.start, gene_location.end)] = {
                    "Gene Name": gene_name,
                    "Gene Length": len(gene_sequence),
                    "location": gene_location,  # Include the gene location in the dictionary
                }
    # Process the pileup file and extract gene information
    mapped_genes_output = []
    with open(pileup_out, "r") as pileup_file:
        for line in pileup_file:
            position, depth = line.strip().split("\t")
            position = int(position)
            depth = int(depth)

            for gene_pos, gene_info in mapped_genes.items():
                gene_start, gene_end = gene_pos

                if gene_start <= position < gene_end:
                    mapped_genes_output.append({
                        "Gene Name": gene_info["Gene Name"],
                        "Gene Position": position,
                        "Depth": depth
                    })
                    break

    # Write the gene information to the CSV file
    csv_output_file = "mapped_genes_info.csv"
    with open(csv_output_file, "w", newline="") as csvfile:
        fieldnames = ["Gene Name", "Gene Position", "Depth"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(mapped_genes_output)

## Merging mapped regions and gene info  
def merge_mapped_and_gene_info(mapped_positions_file, gene_info_file, output_file):
    # Read mapped positions from the BED file
    mapped_regions = []
    with open(mapped_positions_file, "r",newline="") as mapped_file:
        reader = csv.reader(mapped_file)
        next(reader)   ## Skip the header row
        for row in reader: 
            read_name,start_pos,end_pos = row
            mapped_regions.append((read_name, int(start_pos), int(end_pos)))

    # Read gene information from the CSV file generated in extracted_gene_from_mapped 
    gene_info_list= []
    with open(gene_info_file,"r", newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row["Gene Name"]
            gene_pos = int(row["Gene Position"])   ## Convert to integer
            depth = int(row["Depth"])
            gene_info_list.append((gene_name,gene_pos,depth))

# # Create a dictionary to store lowest and highest positions for each gene combination
    gene_positions = {}

    # Iterate through the mapped regions and update the lowest and highest positions
    for read_name, start_pos, end_pos in mapped_regions:
        matched_genes = set()
        highest_depth = 0
        for gene_name, gene_pos, depth in gene_info_list:
            if start_pos <= gene_pos <= end_pos:
                matched_genes.add(gene_name)
                highest_depth = max(highest_depth, depth)

        if matched_genes:
            key = f"{','.join(matched_genes)}_{highest_depth}"
            if key not in gene_positions:
                gene_positions[key] = {"start_pos": start_pos, "end_pos": end_pos}

            # Update the lowest and highest positions if necessary
            gene_positions[key]["start_pos"] = min(gene_positions[key]["start_pos"], start_pos)
            gene_positions[key]["end_pos"] = max(gene_positions[key]["end_pos"], end_pos)

    # Create lists to store the final output data
    start_positions = []
    end_positions = []
    gene_names = []
    depths = []

    # Populate the final output lists
    for key, value in gene_positions.items():
        gene_names.append(','.join(key.split('_')[:-1]))  # Extract gene names
        depths.append(int(key.split('_')[-1]))  # Extract depth
        start_positions.append(value["start_pos"])
        end_positions.append(value["end_pos"])

    # Create the final output DataFrame
    final_output_df = pd.DataFrame({
        'Start_Position': start_positions,
        'End_Position': end_positions,
        'Mapped_Region': [f"{start}-{end}" for start, end in zip(start_positions, end_positions)],
        'Gene_Names': gene_names,
        'Depth': depths
    })

    # Sort the DataFrame by the 'Depth' column in descending order
    final_output_df = final_output_df.sort_values(by='Depth', ascending=False)

    # Save the final output DataFrame to the output file
    final_output_df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(description="Extract Non-Targeted reads from reference FASTA file.")
    parser.add_argument("-targ", help="Path to the input targeted_header.csv file")
    parser.add_argument("-sam", help="Path to the sample.fasta / PC.fasta file")
    parser.add_argument("-ref", help="Path to the reference sequence file")
    parser.add_argument("-gb", help="Path to the genbank file")
    parser.add_argument("-targ", help="Path to the output targeted.txt file")
    parser.add_argument("-non", help="Path to the output non-targeted FASTA file")
    parser.add_argument("-fin", help="Path to the final mapped regions output csv")

    args = parser.parse_args()
    ## Reference for Mapping and Annotation 
    seq = args.ref # Path to the ref genome to which u want to map the sequence  
    ref_genbank_file = args.gb  # Replace with the path to your GenBank file
    # If no arguments are provided, display the help message and exit.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    pathA = os.path.split(args.non)[0]
    print("output folder path: ", pathA)
    
    # Step 1: Read ARG headers from input.tsv and 16S headers from input.csv
    arg_headers = read_csv_file(args.targ)
    targeted_headers = list(set(arg_headers ))
    print("Got the targeted headers")

    # Step 2: Compare with ref.fasta and extract non-matching headers and sequences
    ref_headers, ref_sequences = compare_with_ref(args.sam, targeted_headers)
    non_targeted_headers = [header for header in ref_headers if header not in targeted_headers]
    print("Non_targeted.txt created ! ")

    # Step 3: Write non-matching headers and sequences to non_targeted_read.fasta
    with open(args.non, "w") as file:
        for header, sequence in zip(ref_headers, ref_sequences):
            if header in non_targeted_headers:
                file.write(f">{header}\n")
                file.write(f"{sequence}\n")
    
    # Step 4: Index the ref genome Map Non-Targeted sequence using Minimap2
    ## Outputs 
    index_out_file = os.path.join(pathA,"ref.mmi")
    sam_out_file=os.path.join(pathA,"non_targeted.sam") 
    # Index the reference genome 
    minimap_index_command = f"minimap2 -x map-ont -d {index_out_file} {seq}"
    subprocess.run(minimap_index_command, shell=True, check=True)
    print("Reference genome Indexing done")
    # Map non-targeted sequence to the reference genome 
    minimap_map_command = f"minimap2 -ax map-ont {seq} {args.non} > {sam_out_file}"
    subprocess.run(minimap_map_command, shell=True, check=True)
    print("Mapping Done")

    # Step 5: using Samtools 
    # Convert SAM to BAM 
    bam_out_file=os.path.join(pathA,"non_targeted.bam")
    samtool_command = f"samtools view -S -b {sam_out_file} >{bam_out_file}"
    subprocess.run(samtool_command,shell=True, check =True)
    print("SAM to BAM conversion done")
    #Sort BAM 
    sorted_bam_out=os.path.join(pathA,"non_targeted_sorted.bam")
    samtool_sort_command = f"samtools sort {bam_out_file} -o {sorted_bam_out}"
    subprocess.run(samtool_sort_command,shell =True, check=True)
    print("BAM Sorting DONE !!")
    #Index sorted BAM 
    samtools_index_command = f"samtools index {sorted_bam_out}"
    subprocess.run(samtools_index_command, shell=True, check=True)
    print("BAM_sort index DONE !!")
    #Create pileup file
    pileup_out = os.path.join(pathA,"non_targeted_count.tsv")
    pileup_command = f"samtools mpileup -f {seq} {sorted_bam_out} | awk 'BEGIN {{FS=\"\t\"}} {{print $2\"\t\"$4}}' > {pileup_out}"
    subprocess.run(pileup_command, shell=True, check=True)
    print ("pileup file created")

    # Step 6: Convert BAM to BED 
    BED_out = os.path.join(pathA,"non_targeted_sorted.bed")
    Bed_command = f"bedtools bamtobed -i {sorted_bam_out} > {BED_out}"
    subprocess.run(Bed_command,shell=True, check=True)
    print("BAM to BED done !!")

    # Step 7: Extract mapped positions from BED file
    extract_mapped_regions(BED_out)
    
    # Step 8: Extract gene information and generate CSV file
    extract_gene_from_mapped(ref_genbank_file, pileup_out)
    print("Extracted_Gene name and Depth !! ")

    # Step 9: Merge mapped positions and gene information
    merge_mapped_and_gene_info("mapped_positions.csv", "mapped_genes_info.csv",args.fin)

    print("Process Completed!")

if __name__ == "__main__":
    main()

