
#!/bin/bash

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_fasta_file"
    exit 1
fi

# Define input filename
input_fasta=$1
input_filename=$(basename "$input_fasta" .fasta)

# Define output directory
output_dir="${input_filename}_output"

# Check if output directory already exists
if [ -d "$output_dir" ]; then
    echo "Output directory '$output_dir' already exists. Deleting..."
    rm -rf "$output_dir"
    echo "Existing output directory deleted."
fi

mkdir -p "${output_dir}"

# Define log file
log_file="${output_dir}/pipeline_log.txt"
exec > "${log_file}" 2>&1

# Define intermediate filenames
blast_output="${output_dir}/blast_output.txt"       # Define the path for the BLAST output file
filtered_blast_output="${output_dir}/filtered_blast_output.txt"  # Define the path for the filtered BLAST output file
cd_hit_output="${output_dir}/cd_hit_output.fasta"   # Define the path for the CD-HIT output file

# Ensure the directories exist
mkdir -p "$(dirname "${blast_output}")"
mkdir -p "$(dirname "${filtered_blast_output}")"
mkdir -p "$(dirname "${cd_hit_output}")"

# Function to check if a file is empty
check_file_empty() {
    if [ ! -s "$1" ]; then
        echo "Error: $1 is empty. Exiting."
        exit 1
    fi
}

# Run all-against-all BLAST
echo "Running all_vs_all blast..."
blastp -num_threads 20 -query "${input_fasta}" -subject "${input_fasta}" -evalue 1e-10 -comp_based_stats F -ungapped -max_target_seqs 500 -threshold 11 -outfmt "6 sseqid length evalue sseq sstart send qcovhsp" -out "${blast_output}" || exit 1
check_file_empty "${blast_output}"

# Filter BLAST output for sequences with coverage less than 50%
awk '$7 < 30' "${blast_output}" > "${filtered_blast_output}" || exit 1
check_file_empty "${filtered_blast_output}"

# Extract sequences from the filtered BLAST output
awk '{if($2>10 && $2<500){print ">" $1 "/" $5 "~" $6 "\n" $4}}' "${filtered_blast_output}" | sed 's/-//g' > "${cd_hit_output}" || exit 1
check_file_empty "${cd_hit_output}"

# Cluster sequences
echo "Clustering fragments..."
cd-hit  -i "${cd_hit_output}" -o "${output_dir}/${input_filename}_95_clustering" -c 0.99 -aL 0.6 -aS 0.95 -n 2 -T 20 || exit 1
mmseqs easy-cluster "${output_dir}/${input_filename}_95_clustering" "${output_dir}/${input_filename}_0.2_clust_DB" tmp -c 0.9 --cov-mode 3 --min-seq-id 0.3 -s 7.5 --cluster-reassign 1 --cluster-steps 5 --threads 20 || exit 1

# Generate list of sequences with counts
echo "Generating list of sequences with counts..."
awk '{print $0}' "${output_dir}/${input_filename}_0.2_clust_DB_all_seqs.fasta" | uniq -c > "${output_dir}/list_of_all_with_counts" || exit 1

# Display the content of list_of_all_with_counts
echo "Contents of list_of_all_with_counts:"
cat "${output_dir}/list_of_all_with_counts"

# Create clusters directory
clusters_dir="${output_dir}/${input_filename}_2.0_clusters"
mkdir -p "${clusters_dir}"

# Split the list of sequences into clusters
echo "Splitting sequences into clusters..."
csplit --digits=3 --quiet --prefix="${clusters_dir}/domain" "${output_dir}/list_of_all_with_counts" "/2 >/" "{*}" || exit 1
sed -i 's/      2 //' "${clusters_dir}/domain"* || exit 1
sed -i 's/      1 //' "${clusters_dir}/domain"* || exit 1

# List files in clusters directory
echo "Files in ${clusters_dir}:"
ls -lh "${clusters_dir}"

# Remove small clusters
echo "Removing small clusters..."
for file in "${clusters_dir}/domain"*; do
    NUM=$(grep -c ">" "$file")
    if [ $NUM -lt 5 ]; then
        rm -f "$file"
    fi
done

# Compile HMMs for each cluster
echo "Compiling HMMs..."
MSAs_dir="${output_dir}/MSAs"
HMMs_dir="${output_dir}/HMMs"
mkdir -p "${MSAs_dir}"
mkdir -p "${HMMs_dir}"

# Only process files that match the domain pattern to avoid processing unintended files
for file in "${clusters_dir}/domain"*; do
    if [ -f "$file" ]; then
        echo "Running clustalo on $file"
        clustalo -i "$file" -o "${MSAs_dir}/$(basename "$file")" --threads=20 || exit 1
    fi
done

# List files in MSAs_dir
echo "Files in ${MSAs_dir}:"
ls -lh "${MSAs_dir}"

for file in "${MSAs_dir}/domain"*; do
    if [ -f "$file" ]; then
        echo "Running hmmbuild on $file"
        hmmbuild --cpu 20 "${HMMs_dir}/$(basename "$file" .fasta).hmm" "$file" || exit 1
    fi
done

# List files in HMMs_dir
echo "Files in ${HMMs_dir}:"
ls -lh "${HMMs_dir}"

# Check if HMM files are created before creating the HMM database
if [ -z "$(ls -A ${HMMs_dir})" ]; then
    echo "Error: No HMM files created. Exiting."
    exit 1
fi

# Create HMM database
echo "Creating HMM database..."
cat "${HMMs_dir}"/*.hmm > "${output_dir}/${input_filename}_clusters_HMM_DB" || exit 1
hmmpress "${output_dir}/${input_filename}_clusters_HMM_DB" || exit 1

echo "Running HMMscan against custom HMM database..."
hmmscan --cpu 20 --domT 50 --domtblout "${output_dir}/${input_filename}_DB_hmmscan_out" "${output_dir}/${input_filename}_clusters_HMM_DB" "${input_fasta}" || exit 1

# Run Python script to reduce overlapping domains
python domain_grouping.py "${output_dir}/${input_filename}_DB_hmmscan_out" "${output_dir}/selected_domains.txt"
echo "Number of selected domains: $(wc -l < "${output_dir}/selected_domains.txt")"

# Create a new directory for selected domains' HMMs
selected_hmms_dir="${output_dir}/selected_hmms"
mkdir -p "${selected_hmms_dir}"

# Transfer selected domains' HMMs to the new directory
for domain in $(cat "${output_dir}/selected_domains.txt"); do
    cp "${HMMs_dir}/${domain}.hmm" "${selected_hmms_dir}/"
done

# Refine HMMs with HMMalign and rebuild HMMs
echo "Refining HMMs with HMMalign and rebuilding..."
for hmm_file in "${selected_hmms_dir}"/*.hmm; do
    hmm_name=$(basename "$hmm_file" .hmm)
    alignment_file="${output_dir}/HMMalign/${hmm_name}.sto"
    hmmbuild_output_file="${selected_hmms_dir}/${hmm_name}.hmm"
    
    mkdir -p "${output_dir}/HMMalign"
   
    echo "Running hmmalign on $hmm_file"
    hmmalign --outformat stockholm --trim "$hmm_file" "$input_fasta" > "$alignment_file" || exit 1
    
    echo "Rebuilding HMM with alignment from $alignment_file"
    hmmbuild --fragthresh 0.5 "$hmmbuild_output_file" "$alignment_file" || exit 1
done

# Concatenate selected domains' HMMs into a new HMM database
cat "${selected_hmms_dir}"/*.hmm > "${output_dir}/${input_filename}_selected_clusters_HMM_DB" || exit 1
hmmpress "${output_dir}/${input_filename}_selected_clusters_HMM_DB" || exit 1

# Run HMMscan against the selected domains' HMM database
echo "Running HMMscan against selected domains' HMM database..."
hmmscan --cpu 20 -o "${output_dir}/${input_filename}_selected_DB_hmmscan_out" "${output_dir}/${input_filename}_selected_clusters_HMM_DB" "${input_fasta}" || exit 1

# Run HMMscan against Pfam database
echo "Running HMMscan against Pfam database..."
hmmscan --cpu 20 -o "${output_dir}/${input_filename}_pfam_hmmscan_out" /sci/labs/asafle/alexlevylab/icore-data/new_pfam_db/Pfam-A.hmm "${input_fasta}" || exit 1

# Run CATH on Pfam results
echo "Running CATH on Pfam results..."
cath-resolve-hits.ubuntu-20.04 --input-format hmmscan_out --worst-permissible-bitscore 10 --hits-text-to-file "${output_dir}/${input_filename}_pfam_cath_resolved_hits" --html-exclude-rejected-hits --html-output-to-file "${output_dir}/${input_filename}_pfam_hits_graphic.html" "${output_dir}/${input_filename}_pfam_hmmscan_out"

# Run CATH on selected domains results
echo "Running CATH on selected domains results..."
cath-resolve-hits.ubuntu-20.04 --input-format hmmscan_out --worst-permissible-bitscore 25 --hits-text-to-file "${output_dir}/${input_filename}_selected_cath_resolved_hits" --html-exclude-rejected-hits --html-output-to-file "${output_dir}/${input_filename}_selected_hits_graphic.html" "${output_dir}/${input_filename}_selected_DB_hmmscan_out"

# Rank resolved domains by prevalence in hit pool
echo "Ranking resolved domains by prevalence..."
awk -F " " '{print $2}' "${output_dir}/${input_filename}_cath_resolved_hits" | sort | uniq -c | sort -rn > "${output_dir}/${input_filename}_prevalent_hits"

# Run the selected_domains_extraction.sh script
echo "Running selected domains extraction, consensus generation, and tripling..."
./selected_domains_extraction.sh "${MSAs_dir}" "${output_dir}/selected_domains.txt" 2

# Run HMMscan with selected domains database against cleaned consensus lines
echo "Running HMMscan with selected domains database against cleaned consensus lines..."
hmmscan --domtblout "${output_dir}/cleaned_consensus_hmmscan.domtblout" \
        --cpu 20 \
        --domT 120 \
        --noali \
        "${output_dir}/${input_filename}_selected_clusters_HMM_DB" \
        "${output_dir}/cleaned_sequences.fasta" \
        > "${output_dir}/cleaned_consensus_hmmscan.out" || exit 1

# Check if the domtblout file was created successfully
if [ ! -f "${output_dir}/cleaned_consensus_hmmscan.domtblout" ]; then
    echo "Error: HMMscan failed to produce the domtblout file."
    exit 1
fi

echo "HMMscan completed. Results saved in ${output_dir}/cleaned_consensus_hmmscan.domtblout"

# parsing script to determine HMM groups
echo "Parsing HMMscan results..."
python parse_hmmscan_results.py "${output_dir}/cleaned_consensus_hmmscan.domtblout" "${output_dir}/parsed_hmmscan_results.txt"

# Define paths and parameters
DOMAIN_FILES="${output_dir}/${input_filename}_selected_cath_resolved_hits,${output_dir}/${input_filename}_pfam_cath_resolved_hits"
FASTA_FILE="${input_fasta}"
OUTPUT_FILE="${output_dir}/${input_filename}_domain_plots.html"
SCORE_THRESHOLD=10

# Run the domain plotting script
python domain_plot_try1.py "${DOMAIN_FILES}" "${FASTA_FILE}" "${OUTPUT_FILE}" --score_threshold ${SCORE_THRESHOLD}

echo "Domain plots generated and saved to ${OUTPUT_FILE}"


# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "${output_dir}/${input_filename}_95_clustering"
rm -rf tmp

echo "Pipeline completed successfully."
