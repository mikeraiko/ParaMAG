rule all:
    input:
        "qc_reports/fastqc_done.txt",
        "spades_output/contigs.fasta",
        "maxbin_output/contig_bins.fasta",
        "blobtools_output/plots/blobtools_plot.png",
        "selected_bin/fungal_bin.fasta",
        "busco_output/fungal/fungi_full_table.tsv",
        "busco_output/parasitic_fungi/parasitic_fungi_full_table.tsv",
        "super_tree_output.tre"

rule fastqc:
    input:
        "input_reads/sample_R1.fastq", "input_reads/sample_R2.fastq"
    output:
        "qc_reports/fastqc_done.txt"
    shell:
        """
        fastqc {input} -o qc_reports/
        touch {output}
        """

rule trimmomatic:
    input:
        r1="input_reads/sample_R1.fastq",
        r2="input_reads/sample_R2.fastq"
    output:
        r1_paired="trimmed/sample_R1_paired.fastq",
        r1_unpaired="trimmed/sample_R1_unpaired.fastq",
        r2_paired="trimmed/sample_R2_paired.fastq",
        r2_unpaired="trimmed/sample_R2_unpaired.fastq"
    params:
        adapter_file="adapters.fa"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule spades:
    input:
        r1="trimmed/sample_R1_paired.fastq",
        r2="trimmed/sample_R2_paired.fastq"
    output:
        "spades_output/contigs.fasta"
    shell:
        "spades.py --sc -1 {input.r1} -2 {input.r2} -o spades_output/"

rule maxbin:
    input:
        contigs="spades_output/contigs.fasta",
        abund="abundance_data.tsv"
    output:
        "maxbin_output/contig_bins.fasta"
    shell:
        "run_MaxBin.pl -contig {input.contigs} -out maxbin_output -abund {input.abund}"

rule blobtools:
    input:
        contigs="maxbin_output/contig_bins.fasta",
        coverage="maxbin_output/coverage.tsv",
        taxonomy="taxonomy.tsv"
    output:
        "blobtools_output/plots/blobtools_plot.png"
    shell:
        """
        blobtools create -i {input.contigs} -c {input.coverage} -t {input.taxonomy} --taxrule best --out blobtools_output/
        blobtools plot -i blobtools_output.blobDB.json -o blobtools_output/plots/
        """

rule select_fungal_bin:
    input:
        "blobtools_output/blobtools_plot.png"
    output:
        "selected_bin/fungal_bin.fasta"
    shell:
        """
        # Extract the bin corresponding to fungi using BlobTools filtering based on taxonomic classification.
        blobtools filter -i blobtools_output.blobDB.json -x Fungi -o selected_bin/fungal_bin.fasta
        """

rule busco_fungi:
    input:
        "selected_bin/fungal_bin.fasta"
    output:
        "busco_output/fungal/fungi_full_table.tsv"
    params:
        lineage="fungi_odb10"
    shell:
        "busco -i {input} -l {params.lineage} -o busco_output/fungal -m genome"

rule busco_parasitic_fungi:
    input:
        "selected_bin/fungal_bin.fasta"
    output:
        "busco_output/parasitic_fungi/parasitic_fungi_full_table.tsv"
    params:
        lineage="parasitic_fungi_odb10"  # Assumes a custom lineage for parasitic fungi is available.
    shell:
        "busco -i {input} -l {params.lineage} -o busco_output/parasitic_fungi -m genome"

rule phylogenomics:
    input:
        "busco_output/fungal/fungi_full_table.tsv"
    output:
        "super_tree_output.tre"
    shell:
        """
        mafft --auto selected_bin/fungal_bin.fasta > aligned_sequences.fasta
        raxmlHPC -s aligned_sequences.fasta -n phylogenetic_tree -m PROTGAMMAJTT -p 12345
        java -jar astral.jar -i gene_trees.tre -o super_tree_output.tre
        """

