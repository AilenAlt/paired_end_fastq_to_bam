configfile: 'config.yaml'

SRRs = [
    "SRR10485675",
    "SRR10485677",
    "SRR10485678",
]

rule all:
    input:
        expand('align/{srr}_filtered.bam', srr=SRRs)

rule trim:
    input: 
        in1 = '{srr}_1.fastq.gz',
        in2='{srr}_2.fastq.gz'
    output: 
        out = '{srr}_1_trimm.fastq.gz',
        out2 = '{srr}_2_trimm.fastq.gz'
    params:
        adapters = config['Adapters'],
        out = '{srr}_1_trimm.fastq',
        out2 = '{srr}_2_trimm.fastq'
    log:
        '{srr}_trimm.log'
    shell:
        """
        bbduk.sh in={input.in1} in2={input.in2} ref={params.adapters} ktrim=r k=23 hdist=1 tbo tpe mink=11 mlf=.5 rcomp=t out={params.out} out2={params.out2} 2> {log} && gzip {params.out} {paramsS.out2}
        """

rule qc:
    input:
        in1 = '{srr}_1_trimm.fastq.gz',
        in2='{srr}_2_trimm.fastq.gz'
    output:
        out = 'qc/{srr}_1_trimm_fastqc.html', 
        out2 = 'qc/{srr}_2_trimm_fastqc.html' 
    threads:
        8
    shell:
        """
        fastqc -o ./qc --threads {threads} {input.in1} {input.in2}
        """

rule align:
    input:
        fastq1 = '{srr}_1_trimm.fastq.gz',
        qc1 = 'qc/{srr}_1_trimm_fastqc.html',
        fastq2 = '{srr}_2_trimm.fastq.gz',
        qc2 = 'qc/{srr}_2_trimm_fastqc.html'
    output:
        'align/{srr}.sam'
    threads:
        8
    conda:
        'bowtie'
    log:
        '{srr}_bowtie.log'
    shell:
        """
        bowtie2 -x dm6_index -1 {input.fastq1} -2 {input.fastq2} --threads {threads} -S {output} 2> {log} && rm {input.fastq1} {input.fastq2}
        """

rule sam_to_bam:
    input:
        'align/{srr}.sam'
    output:
        'align/{srr}_sorted.bam'
    threads:
        8
    conda:
        'samtools'
    shell:
        """
        samtools view -b --threads {threads} {input} | samtools sort --threads {threads} -o {output} - && rm {input}
        """

rule mark_duplicates:
    input:
        'align/{srr}_sorted.bam'
    output:
        'align/{srr}_flagged.bam'
    conda:
        'picard'
    log:
        '{srr}_duplicates.log'
    shell:
        """
        picard MarkDuplicates -I {input} -M {log} -O {output} && rm {input}
        """

rule filter_alignment:
    input:
        'align/{srr}_flagged.bam'
    output:
        'align/{srr}_filtered.bam'
    conda:
        'samtools'
    threads:
        8
    shell:
        """
        samtools view -hb {input} --threads {threads} -q 30 -F 3332 -f 3 -o {output}
        """
