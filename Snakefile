SAMPLES = ['SPR1-1_bulkmetagenome_S42', 'SPR1-1_virome_S35',
           'SPR2-1_bulkmetagenome_S34', 'SPR2-1_virome_S46',
           'SPR3-3_bulkmetagenome_S48', 'SPR3-3_virome_S37',
           'SPR4-2_bulkmetagenome_S53', 'SPR4-2_virome_S27',
           'SPR5-3_bulkmetagenome_S52', 'SPR5-3_virome_S55']

VIRS = ["EarthsVirome_31697", "EarthsVirome_41546", "EarthsVirome_48981"]

rule all:
    input:
#        expand("aggregated_checkpoints/aggregate_sgc_alignments/{sample}.txt", sample = SAMPLES),
        expand('outputs/multiqc_flagstat/{sample}_multiqc.html', sample = SAMPLES),
        expand("aggregated_checkpoints/aggregate_sgc_assembly/{sample}.txt", sample = SAMPLES)

rule download_adapters:
    output: "inputs/adapters.fa"
    shell:'''
    wget -O {output} https://raw.githubusercontent.com/dib-lab/dib_rotation/master/_static/adapters.fa
    '''

rule adapter_trim_reads:
    input:
        r1 = "inputs/annie/{sample}_L001_R1_001.fastq.gz", 
        r2 = "inputs/annie/{sample}_L001_R2_001.fastq.gz",
        adapters = 'inputs/adapters.fa'
    output:
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{sample}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{sample}_o2.trim.fq.gz'
    conda: 'trim.yml'
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:31  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''


rule kmer_trim_reads:
    input: 
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'trim.yml'
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 30e9 -V - -o {output}
    '''

rule split_fasta_to_contigs:
    output: "outputs/contigs/{virs}.fa"
    input: "inputs/queries/recovered_own_data.fa"
    conda: 'fasplit.yml'
    shell: '''
    faSplit byname {input} {output}
    '''

rule index_virs:
    output: "outputs/contigs/{virs}.fa.bwt"
    input: "outputs/contigs/{virs}.fa"
    conda: 'bwa.yml'
    shell:'''
    bwa index {input}
    '''

checkpoint sgc_vir_queries:
    input: 
        conf = "inputs/conf/{sample}_conf.yml",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        virs = expand("outputs/contigs/{virs}.fa", virs = VIRS) 
    conda:  "spacegraphcats.yml"
    output: directory("outputs/sgc/{sample}_k31_r1_search_oh0")
    params: outdir = "outputs/sgc"
    shell:'''
    python -m spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} 
    '''

rule align_reads_to_queries:
    input: 
        sgc_reads = "outputs/sgc/{sample}_k31_r1_search_oh0/{genome}.fa.cdbg_ids.reads.fa.gz",
        index = 'outputs/contigs/{genome}.fa.bwt',
        virs = 'outputs/contigs/{genome}.fa'
    output: "outputs/sgc_align/{sample}/{genome}.sam"
    conda: "bwa.yml"
    shell:'''
    bwa mem {input.virs} {input.sgc_reads} > {output}
    '''

rule flagstat:
    input: 'outputs/sgc_align/{sample}/{genome}.sam'
    output: 'outputs/flagstat/{sample}/{genome}.flagstat'
    conda: "samtools.yml"
    shell:'''
    samtools flagstat {input} > {output}
    '''

def aggregate_sgc_vir_queries(wildcards):
    checkpoint_output = checkpoints.sgc_vir_queries.get(**wildcards).output[0]    
    file_names = expand("outputs/flagstat/{sample}/{genome}.flagstat",
                        sample = wildcards.sample, 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.fa.cdbg_ids.reads.fa.gz")).genome)
    return file_names

rule multiqc_flagstat:
    input: aggregate_sgc_vir_queries
    #'outputs/flagstat/{sample}/{genome}.flagstat'
    output: 'outputs/multiqc_flagstat/{sample}_multiqc.html'
    params: 
        outdir = "outputs/multiqc_flagstat",
        indir = lambda wildcards: "outputs/flagstat/" + wildcards.sample
    conda: "multiqc.yml"
    shell:'''
    multiqc {params.indir} -o {params.outdir} -n {wildcards.sample}_multiqc
    '''

#rule aggregate_alignments:
#    input: aggregate_sgc_vir_queries
#    output: "aggregated_checkpoints/aggregate_sgc_alignments/{sample}.txt"
#    shell:'''
#    touch {output}
#    '''

rule sam_unmapped_reads:
    input: "outputs/sgc_align/{sample}/{genome}.sam"
    output: "outputs/sgc_align/{sample}/{genome}_unmapped.sam"
    conda: 'samtools.yml'
    shell:'''
    samtools view -f 4 {input} > {output}
    ''' 

rule unmapped_reads_to_fastq:
    input: "outputs/sgc_align/{sample}/{genome}_unmapped.sam"
    output: "outputs/sgc_align/{sample}/{genome}_unmapped.fa"
    conda: 'samtools.yml'
    shell:'''
    samtools fasta {input} > {output}
    '''

rule megahit_unmapped_reads:
    input: "outputs/sgc_align/{sample}/{genome}_unmapped.fa"
    output: "outputs/megahit_unmapped/{sample}/{genome}.contigs.fa"
    conda: 'megahit.yml'
    shell:'''
    megahit -r {input} --min-contig-len 142 \
        --out-dir {wildcards.sample}_{wildcards.genome}_megahit \
        --out-prefix {wildcards.sample}_{wildcards.genome}
    mv {wildcards.sample}_{wildcards.genome}_megahit/{wildcards.sample}_{wildcards.genome}.contigs.fa {output}
    rm -rf {wildcards.sample}_{wildcards.genome}_megahit
    '''

def aggregate_sgc_vir_queries2(wildcards):
    checkpoint_output = checkpoints.sgc_vir_queries.get(**wildcards).output[0]    
    file_names = expand("outputs/megahit_unmapped/{sample}/{genome}.contigs.fa",
                        sample = wildcards.sample, 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.fa.cdbg_ids.reads.fa.gz")).genome)
    return file_names

rule aggregate_assembly:
    input: aggregate_sgc_vir_queries2
    output: "aggregated_checkpoints/aggregate_sgc_assembly/{sample}.txt"
    shell:'''
    touch {output}
    '''
