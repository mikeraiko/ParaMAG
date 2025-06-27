configfile: "config.yaml"

reads1 = config["reads1"]
reads2 = config["reads2"]
adapter_file = config["adapter_file"]

rule all:
    input:
        "qc_reports/fastqc_done.txt",
        "spades_output/contigs.fasta",
        "abundance_data.tsv",
        "maxbin_output/contig_bins.fasta",
        "coverage.tsv",
        "blobtools_output/plots/blobtools_plot.png",
        "selected_bin/fungal_bin.fasta",
        "busco_output/fungal/fungi_full_table.tsv",
        "busco_output/parasitic_fungi/parasitic_fungi_full_table.tsv",
        "super_tree_output.tre",
        "report/summary_report.md"

rule fastqc:
    input:
        reads1, reads2
    output:
        "qc_reports/fastqc_done.txt"
    shell:
        """
        fastqc {input} -o qc_reports/
        touch {output}
        """

rule trimmomatic:
    input:
        r1=reads1,
        r2=reads2
    output:
        r1_paired="trimmed/sample_R1_paired.fastq",
        r1_unpaired="trimmed/sample_R1_unpaired.fastq",
        r2_paired="trimmed/sample_R2_paired.fastq",
        r2_unpaired="trimmed/sample_R2_unpaired.fastq"
    params:
        adapter_file=adapter_file
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule spades:
    input:
        r1="trimmed/sample_R1_paired.fastq",
        r2="trimmed/sample_R2_paired.fastq"
    output:
        "spades_output/contigs.fasta"
    shell:
        "spades.py --sc -1 {input.r1} -2 {input.r2} -o spades_output/"

rule abundance:
    input:
        contigs="spades_output/contigs.fasta",
        r1="trimmed/sample_R1_paired.fastq",
        r2="trimmed/sample_R2_paired.fastq"
    output:
        "abundance_data.tsv"
    shell:
        """
        bowtie2-build {input.contigs} contigs_index
        bowtie2 -x contigs_index -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o aln.bam
        samtools index aln.bam
        jgi_summarize_bam_contig_depths --outputDepth {output} aln.bam
        """

rule coverage:
    input:
        bam="aln.bam"
    output:
        "coverage.tsv"
    shell:
        "samtools idxstats {input.bam} > {output}"

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
        coverage="coverage.tsv",
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
        lineage="parasitic_fungi_odb10"
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

rule report:
    input:
        busco="busco_output/fungal/fungi_full_table.tsv",
        tree="super_tree_output.tre"
    output:
        "report/summary_report.md"
    shell:
        """
        echo "# Pipeline Summary Report" > {output}
        echo "\n## BUSCO Results" >> {output}
        tail -n 20 {input.busco} >> {output}
        echo "\n## Phylogenetic Tree File" >> {output}
        echo "Tree saved as: {input.tree}" >> {output}
        """
