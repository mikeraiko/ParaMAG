# Snakefile (с использованием IQ-TREE вместо RAxML)
configfile: "config.yaml"

reads1 = config["reads1"]
reads2 = config["reads2"]
adapter_file = config["adapter_file"]

bins = glob_wildcards("maxbin_output/bin.{binid}.fasta")

rule all:
    input:
        "qc_reports/fastqc_done.txt",
        "spades_output/contigs.fasta",
        "abundance_data.tsv",
        expand("maxbin_output/bin.{binid}.fasta", binid=bins.binid),
        expand("busco_output/bin.{binid}/short_summary.specific.fungi_odb10.txt", binid=bins.binid),
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
        bam="aln.bam",
        depth="abundance_data.tsv"
    shell:
        """
        bowtie2-build {input.contigs} contigs_index
        bowtie2 -x contigs_index -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output.bam}
        samtools index {output.bam}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}
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
        expand("maxbin_output/bin.{i:03}.fasta", i=range(1, 6))
    shell:
        """
        run_MaxBin.pl -contig {input.contigs} -out maxbin_output/bin -abund {input.abund}
        """

rule busco:
    input:
        "maxbin_output/bin.{binid}.fasta"
    output:
        "busco_output/bin.{binid}/short_summary.specific.fungi_odb10.txt"
    params:
        lineage="fungi_odb10",
        outdir="busco_output/bin.{binid}"
    shell:
        "busco -i {input} -l {params.lineage} -o {params.outdir} -m genome"

rule align_bins:
    input:
        expand("maxbin_output/bin.{binid}.fasta", binid=bins.binid)
    output:
        "aligned_bins.fasta"
    shell:
        """
        cat {input} > concatenated_bins.fasta
        mafft --auto concatenated_bins.fasta > {output}
        """

rule phylogenomics:
    input:
        "aligned_bins.fasta"
    output:
        "super_tree_output.tre"
    shell:
        """
        iqtree2 -s {input} -m MFP -bb 1000 -nt AUTO -pre iqtree_output
        cp iqtree_output.treefile {output}
        """

rule report:
    input:
        summaries=expand("busco_output/bin.{binid}/short_summary.specific.fungi_odb10.txt", binid=bins.binid),
        tree="super_tree_output.tre"
    output:
        "report/summary_report.md"
    shell:
        """
        echo "# Pipeline Summary Report" > {output}
        echo "\n## BUSCO Summaries" >> {output}
        for f in {input.summaries}; do echo "\n### $$f" >> {output}; tail -n 10 $$f >> {output}; done
        echo "\n## Phylogenetic Tree File" >> {output}
        echo "Tree saved as: {input.tree}" >> {output}
        """

