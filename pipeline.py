#!usr/bin/env python

"""
Title: NGS pipeline for paired reads
Date: 03/03/2022
Author: Jay Miles

Required installs:
    -FastQC        (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    -BWA            (http://bio-bwa.sourceforge.net/)
    -SAMtools       (http://www.htslib.org/)
    -FreeBayes      (https://github.com/freebayes/freebayes)
    -VCFtools       (https://vcftools.github.io/examples.html)
    -BCFtools       (http://www.htslib.org/)
    -ANNOVAR        (https://annovar.openbioinformatics.org/en/latest/)

Required input:
    -Reference genome assembly file
    (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/)
    -A forward reads .fastq file
    -A reverse reads .fastq file
    -.bed file defining genes to call variants in (example: lqts_genes.bed)
    -convert_scaffolsds text file to convert RefSeq chromosome notation to UCSC

Process based on LQTS example: (creates files ~20.7GB, duration ~2h 15m)
    -Generate .fastq file quality metrics using FastQC          (1.1 MB, ~1m)
    -Generate reference genome index using BWA                  (5.7GB, ~1.2h)
    -BWA-MEM alignment of reads to reference genome             (3.7GB, ~30m)
    -Convert .sam file to .bam using SAMtools                   (2.9GB, ~10s)
    -Collate reads using SAMtools                               (2.9GB, ~1m)
    -Add mate scores for markdup using SAMtools                 (910.5MB, ~2m)
    -Sort reads using SAMtools                                  (3.0GB, ~1m)
    -Mark duplicate reads using SAMtools                        (593.9MB, ~2m)
    -Call variants using Freebayes                              (54.6MB, ~15m)
    -Find and replace RefSeq chromosome notations               (53.6MB, ~3s)
    -Filter against QUAL scores and .bed file using VCFtools    (41.2kB, <1s)
    -Sort .vcf file using bcftools                              (41.2kB, <1s)
    -Download ANNOVAR annotation databases                      (1.1GB, 1-3m)
    -Convert .vcf file to ANNOVAR file format                   (30.1kB, <1s)
    -Annotate variants using ANNOVAR                            (42.3 kB, ~5s)

"""


import subprocess
from datetime import datetime as dt


def run_fastqc(file):
    """ Get quality metrics for the supplied file using FastQC

    Args:
        file [string]: .fastq/.sam/.bam file to get quality metrics for

    Command line: fastqc --noextract fastqc_reads/cardiac_R1.fastqsanger

    Outputs:
        cardiac_R1.fastqsanger_fastqc.html (approx. 236.4 kB)
        cardiac_R1.fastqsanger_fastqc.zip (approx. 289.8 kB)
        cardiac_R2.fastqsanger_fastqc.html (approx. 237.7 kB)
        cardiac_R2.fastqsanger_fastqc.zip (approx. 292.1 kB)
    """

    start_time = dt.now()
    print('{} Running FastQC'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'fastqc',
        '--noextract',
        '{}'.format(file),
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('FastQC took {}\n'.format(str(duration)))


def bwa_index(file):
    """
    Use 'bwa index' to construct index files for the supplied file
    N.B. V LONG RUNTIME

    Args:
        file [string]: path to .fasta/.sam/.bam file to index

    Command line:               bwa index ref_genome/hg19_assembly
    Runtime (hg19_assembly):    approx. 73 minutes (710 iterations)

    Output (hg19_assembly) (total approx. 5.7 GB):
        hg19_assembly.amb (approx. 9.4 KB)
        hg19_assembly.ann (approx. 35.0 KB)
        hg19_assembly.bwt (approx. 3.2 GB)
        hg19_assembly.pac (approx. 808.7 MB)
        hg19_assembly.sa (approx. 1.6 GB)
    """

    start_time = dt.now()
    print('{} BWA-MEM genome indexing'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bwa',
        'index',
        file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Indexing took {}\n'.format(str(duration)))


def bwa_mem_paired(ref_genome, for_reads, rev_reads, output_file):
    """
    Use BWA-MEM to align reads to reference genome
    N.B. LONG RUNTIME

    Args:
        ref_genome [string]: reference genome assembly file (indexed)
        for_reads [string]: .fastq file with forward reads
        rev_reads [string]: .fastq file with reverse reads
        output_file [string]: path to output .sam file

    Command line:   bwa mem \
                        ref_genome/hg19_assembly \
                        fastq_reads/reads_for.fastqsanger \
                        fastq_reads/reads_rev.fastqsanger \
                        -o bam_files/aln.sam

    Runtime:        approx. 28 minutes
    Output:         aln.sam (approx. 3.7GB)
    """

    start_time = dt.now()
    print('{} BWA-MEM alignment'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bwa',
        'mem',
        ref_genome,
        for_reads,
        rev_reads,
        '-o',
        output_file
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Alignment took {}\n'.format(str(duration)))

    return output_file


def bowtie_index(ref_genome):
    """
    THIS FUNCTION IS NOT CURRENTLY USED

    Use 'bowtie2 index' to construct index files for the supplied file

    Args:
        file [string]: path to .fasta/.sam/.bam file to index

    Command line:   bowtie2-build \
                        -f --large-index \
                        ref_genome/hg19_assembly \
                        ref_genome/hg19_assembly

    Runtime (hg19_assembly): approx. 69 minutes

    Output (hg19_assembly) (total approx. 748 MB):
        hg19_assembly.1.bt2l (20.6 kB)
        hg19_assembly.2.bt2l (4 B)
        hg19_assembly.3.bt2l (13.3 kB)
        hg19_assembly.4.bt2l (747.9 MB)

    """

    start_time = dt.now()
    print('{} BowTie2 genome indexing'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bowtie2-build',
        '-f',  # input is in fasta format
        '--large-index',
        ref_genome,
        ref_genome,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Indexing took {}\n'.format(str(duration)))


def bowtie_alignment(ref_genome, for_reads, rev_reads, output_file):
    """
    THIS FUNCTION IS NOT CURRENTLY USED

    Use BowTie2 to align reads to reference genome

    Args:
        ref_genome [string]: reference genome assembly file (indexed)
        for_reads [string]: .fastq file with forward reads
        rev_reads [string]: .fastq file with reverse reads
        output_file [string]: path to output .sam file

    Command line:   bowtie2 \
                        -q \
                        -x ref_genome/hg19_assembly \
                        -1 fastq_reads/reads_for.fastqsanger \
                        -2 fastq_reads/reads_rev.fastqsanger \
                        -S bam_files/aln.sam

    $BOWTIE2_INDEXES=/home/jay/projects/ngs_pipeline/ref_genome/

    Runtime:        approx.
    Output:         aln.sam (approx. )
    """

    start_time = dt.now()
    print('{} BowTie2 alignment'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bowtie2',
        '-qx',
        ref_genome,
        '-1',
        for_reads,
        '-2',
        rev_reads,
        '-S',
        output_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Alignment took {}\n'.format(str(duration)))

    return output_file


def sam_to_bam(sam_file, output_file):
    """
    Use 'samtools view' to convert aln.sam to .bam format

    Args:
        sam_file [string]: path to input .sam file
        output_file [string]: path to output .bam file

    Command line:   samtools view \
                        -bhuo bam_files/aln.bam \
                        bam_files/aln.sam

    Output:         aln.bam (approx. 2.9GB GZip), ~10s
    """

    start_time = dt.now()
    print('{} converting .sam to .bam'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'samtools',
        'view',
        '-bhuo',  # output an uncompressed .bam file with a header line
        output_file, # output file
        sam_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('.sam conversion took {}\n'.format(str(duration)))

    return output_file


def collate_reads(bam_file, output_file):
    """
    Use 'samtools collate' to collate reads in aln.bam

    Args:
        bam_file [string]: path to input .bam file
        output_file [string]: path to output .bam file

    Command line:   samtools collate \
                        -uo bam_files/collated.bam \
                        bam_files/aln.bam

    Outputs:        collated.bam (approx. 2.9GB GZip), ~1m
    """

    start_time = dt.now()
    print('{} Collating reads'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'samtools',
        'collate',
        '-u',
        '-o',
        output_file,
        bam_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Read collation took {}\n'.format(str(duration)))

    return output_file


def fixmates(bam_file, output_file):
    """ Run samtools fixmate on collated.bam to enable later use of samtools
    markdup

    Args:
        bam_file [string]: path to input .bam file
        output_file [string]: path to output .bam file

    Command line:   samtools fixmate -mu \
                        bam_files/collated.bam \
                        bam_files/fixmate.bam

    Outputs:        fixmate.bam (approx. 911MB GZip), ~2m
    """

    start_time = dt.now()
    print('{} Mate fixing'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'samtools',
        'fixmate',
        '-mu',  # add mate scores to enable markdup later, output uncompressed
        bam_file,
        output_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Mate fixing took {}\n'.format(str(duration)))

    return output_file


def sort_bam(bam_file, output_file):
    """
    Use 'samtools sort' to sort fixmate.bam

    Args:
        bam_file [string]: name of input .bam file
        output_file [string]: name for output .bam file

    Command line:   samtools sort \
                        -uo bam_files/sorted.bam \
                        bam_files/fixmate.bam

    Outputs:        sorted.bam (approx. 3 GB GZip), ~1m
    """

    start_time = dt.now()
    print('{} Sorting .bam file'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'samtools',
        'sort',
        '-uo',
        output_file,
        bam_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Sorting took {}\n'.format(str(duration)))

    return output_file


def mark_duplicates(bam_file, output_file):
    """
    Use 'samtools markdup' to mark duplicate reads arising from PCR
    artefacts in sorted.bam

    Args:
        bam_file [string]: name of input .bam file
        output_file [string]: name for output .bam file

    Command line:   samtools markdup \
                        bam_files/sorted.bam \
                        bam_files/markdup.bam

    Outputs:        markdup.bam (approx. 594 MB GZip), ~2m
    """

    start_time = dt.now()
    print('{} Marking duplicates'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'samtools',
        'markdup',
        bam_file,
        output_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Marking duplicates took {}\n'.format(str(duration)))

    return output_file


def call_variants(ref_genome, bam_file, output_file):
    """
    Use FreeBayes to call variants in markdup.bam compared to the
    reference genome assembly

    Args:
        ref_genome [string]: reference genome assembly file (indexed)
        bam_file [string]: name of input .bam file
        output_file [string]: name for output .vcf file

    Command line:   freebayes -f \
                        ref_genome/hg19_assembly \
                        bam_files/markdup.bam \
                        > vcf_files/variants.vcf

    Outputs:        variants.vcf (approx. 54.6 MB), ~15m
    """

    start_time = dt.now()
    print('{} Calling variants'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        "freebayes",
        "-f",
        ref_genome,
        bam_file
        ],
        stdout = open(output_file, 'w'))

    end_time = dt.now()
    duration = end_time - start_time
    print('Variant calling took {}\n'.format(str(duration)))

    return output_file


def chromosome_notation(chr_file, file_to_change):
    """
    Convert RefSeq contig notation into UCSC format (needed for filtering)

    Args:
        chr_file [string]: contains chromosome notation info
        sam_file [string]: name of input .sam file

    Output: variants.vcf (approx. 53.6 MB), ~3s
    """

    start_time = dt.now()
    print('{} Fixing chromosome notation'.format(start_time.strftime('%H:%M')))

    # Generate a list of contig name conversions

    chr_names = []

    with open(chr_file, 'r') as reader:
        lines = reader.readlines()

    for line in lines:
        elements = line.split('\t')

        values = [
            elements[0].strip(),  # RefSeq contig name
            elements[1].strip(),  # UCSC contig name
            elements[2].strip(),  # Alternative sequence name
            ]

        chr_names.append(values)

    # Use sed to find and replace strings

    for names in chr_names:

        # Some contigs don't have a UCSC name, so use the alternative one
        if names[1] == 'na':

            subprocess.run([
                'sed',
                '-i',  # change file in-place
                's/{}/{}/g'.format(names[0], names[2]),
                file_to_change,
                ])

        # Otherwise use the UCSC name
        else:

            subprocess.run([
                'sed',
                '-i',
                's/{}/{}/g'.format(names[0], names[1]),
                file_to_change,
                ])


    end_time = dt.now()
    duration = end_time - start_time
    print('Chromosome notation took {}\n'.format(str(duration)))

    return file_to_change


def filter_vcf(vcf_file, bed_file, minQ, output_prefix):
    """  Filter variants in a .vcf file using QUAL values

    Args:
        vcf_file [string]: .vcf file of variants to be filtered
        bed_file [string]: name of bed file to filter regions on
        minQ [string]: minimum value of QUAL for variants to pass
        output_prefix [string]: prefix for output files

    Command line:   vcftools \
                        --vcf vcf_file/variants.vcf \
                        --bed cardiac_bed.bed \
                        --minQ 20 \
                        --recode \
                        --out vcf_file/filtered

    Outputs:      filtered.vcf (approx. 41.2kB), ~0.2s
    """

    start_time = dt.now()
    print('{} Filtering variants'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'vcftools',
        '--vcf',
        vcf_file,
        '--bed',
        bed_file,
        '--minQ',
        minQ,
        '--recode',
        '--out',
        output_prefix,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Variant filtering took {}\n'.format(str(duration)))

    return '{}.recode.vcf'.format(output_prefix)


def sort_variants(vcf_file, output_file):
    """
    Use 'bcftools sort' to sort a .vcf file

    Args:
        vcf_file [string]: name of input .bam file
        output_file [string]: name for output .bam file

    Command line:   bcftools sort \
                        -o vcf_file/sorted.vcf \
                        vcf_file/filtered.recode.vcf

    Outputs:        sorted.vcf (approx. 41.2 kB), <1s
    """

    start_time = dt.now()
    print('{} Sorting .vcf file'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'bcftools',
        'sort',
        '-o',
        output_file,
        vcf_file,
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Sorting took {}\n'.format(str(duration)))

    return output_file


def annovar_downloads():
    """
    Download ANNOVAR databases for gene- and filter-based annotation.
    DBs are downloaded into a 'humandb' folder.

    ClinVar                 clinvar_20210501    2021 05 07    202.7 MB
    RefGene                 refGene             2021 10 19    308.6 MB
    Splicing predictors     dbscsnv11           2015 12 18    576.2 MB

    Running for just the first 3 takes ~1.5m

    The following 3 DBs are very (impractically) big...

    gnomAD                  gnomad211_genome    2017 09 29    14.2 GB
    dbSNP                   avsnp150            2019 03 23    33.1 GB
    Missense predictors     dbnsfp42c           2021 07 10    41.3 GB

    Command line (insert appropriate db name):
    annotate_variation.pl -downdb -buildver hg19 -webfrom annovar <DB> humandb/
    """

    start_time = dt.now()
    print('{} Downloading ANNOVAR dbs'.format(start_time.strftime('%H:%M')))

    databases = [
        'clinvar_20210501',
        'refGene',
        'dbscsnv11',
        # 'avsnp150',
        # 'gnomad211_genome',
        # 'dbnsfp42c',
        ]

    for db in databases:

        subprocess.run([
            'annotate_variation.pl',
            '-downdb',
            '-buildver',
            'hg19',
            '-webfrom',
            'annovar',
            db,
            'humandb/',
            ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Database download took {}\n'.format(str(duration)))


def annovar_format(vcf_file, output_file):
    """
    Convert a .vcf file to ANNOVAR file format

    Args:
        vcf_file [string]: .vcf file to convert
        output_file [string]: converted file in ANNOVAR file format

    Command line:   convert2annovar.pl \
                        -format vcf4 vcf_file/filtered.recode.vcf \
                        -outfile vcf_file/pre_annovar.vcf \
                        -includeinfo \
                        -withzyg

    Outputs:
        pre_annovar.vcf (approx. 30.1 kB), <1s
    """

    start_time = dt.now()
    print('{} Converting to .avinput'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'convert2annovar.pl',
        '-format',
        'vcf4',
        vcf_file,
        '-outfile',
        output_file,
        '-includeinfo',
        '-withzyg',
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('ANNOVAR conversion took {}\n'.format(str(duration)))

    return output_file


def annotate_variants(vcf_file, output_prefix):
    """  Use ANNOVAR to perform gene- and filter-based variant annotation

    Args:
        vcf_file [string]: .vcf file of variants to be annotated
        output_prefix [string]: prefix for output .vcf file

    Command line:

        table_annovar.pl \
            vcf_file/pre_annotation.vcf \
            humandb/ \
            -buildver hg19 \
            -out vcf_file/post_annotation \
            -remove \
            -protocol refGene,clinvar_20210501,dbscsnv11 \
            -operation g,f,f \
            -nastring .

    Outputs: post_annotation.vcf (approx. 42.3 kB), ~5s
    """

    start_time = dt.now()
    print('{} Annotating variants'.format(start_time.strftime('%H:%M')))

    subprocess.run([
        'table_annovar.pl',
        vcf_file,
        'humandb/',
        '-buildver',
        'hg19',
        '-out',
        output_prefix,
        '-remove',
        '-protocol',
        'refGene,clinvar_20210501,dbscsnv11',
        '-operation',
        'g,f,f',
        '-nastring',
        '.',
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('Variant annotation took {}\n'.format(str(duration)))


def compare_vcfs(vcf_1, vcf_2, output_prefix):
    """  Compare two vcf files to identify common and unique variants

    Args:
        vcf_1 [string]: name of first vcf file
        vcf_2 [string]: name of second vcf file
        output_prefix [string]: prefix for output files

    Command line:
    vcftools --vcf vcf_1 --diff vcf_2 --diff-site --out output_prefix

    Outputs:       (approx. ), ~
    """

    start_time = dt.now()

    print('{} Comparing {} and {}'.format(
        start_time.strftime('%H:%M'),
        vcf_1,
        vcf_2,
        ))

    subprocess.run([
        'vcftools',
        '--vcf',
        vcf_1,  # first input .vcf file
        '--diff',
        vcf_2,  # .vcf file to compare to
        '--diff-site',  # look at common and unique sites
        '--out',
        output_prefix,  # prefix for output files
        ])

    end_time = dt.now()
    duration = end_time - start_time
    print('VCF comparison took {}\n'.format(str(duration)))


def main():

    print('PIPELINE STARTED: {}'.format(dt.now().strftime('%H:%M:%S')))
    pipeline_start = dt.now()


    """ Define inputs as strings of paths to input files """

    ref_genome = 'ref_genome/hg19_assembly'

    for_reads = 'fastq_reads/reads_for.fastqsanger'
    rev_reads = 'fastq_reads/reads_rev.fastqsanger'

    chr_file = 'convert_scaffolds'  # convert RefSeq > UCSC chr notation

    bed_file = 'lqts_genes.bed'
    minQ = '20'  # minimum QUAL score to retain a variant during filtering

    bowtie_vars = 'galaxy_downloads/1_vars_bowtie.vcf'
    bwamem_vars = 'galaxy_downloads/1_vars_bwamem.vcf'
    original_vars = 'galaxy_downloads/1_vars_original_snps.vcf'


    """ Get quality metrics for forward and reverse reads files """

    for reads in [for_reads, rev_reads]:
        run_fastqc(reads)


    """ Create genome assembly index and perform alignment """

    bwa_index(ref_genome)

    sam_file = bwa_mem_paired(
        ref_genome,
        for_reads,
        rev_reads,
        'bam_files/aln.sam')


    """ Convert to .bam format and clean files """

    bam_file = sam_to_bam(sam_file, 'bam_files/1_aln.bam')
    collated_bam = collate_reads(bam_file, 'bam_files/2_collated.bam')
    fixedmates_bam = fixmates(collated_bam, 'bam_files/3_fixmate.bam')
    sorted_bam = sort_bam(fixedmates_bam, 'bam_files/4_sorted.bam')
    markdup_bam = mark_duplicates(sorted_bam, 'bam_files/5_markdup.bam')


    """ Call variants using FreeBayes, filter using VCFTools """

    vcf_file = call_variants(
        ref_genome,
        markdup_bam,
        'vcf_files/1_variants.vcf')

    fixed_chrs = chromosome_notation(chr_file, vcf_file)

    filtered_vcf = filter_vcf(
        fixed_chrs,
        bed_file,
        minQ,
        'vcf_files/2_filtered')

    sorted_vcf = sort_variants(filtered_vcf, 'vcf_files/3_sorted.vcf')


    """ Download annotation databases and annotate variants """

    annovar_downloads()

    pre_annotation_vcf = annovar_format(
        sorted_vcf,
        'vcf_files/4_pre_annotation.vcf')

    annotate_variants(pre_annotation_vcf, 'vcf_files/post_annotation')


    pipeline_end = dt.now()
    pipeline_duration = pipeline_end - pipeline_start
    print('PIPELINE TOOK {}'.format(str(pipeline_duration)))


if __name__ == '__main__':
    main()
