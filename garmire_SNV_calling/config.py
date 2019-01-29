"""
config file for SNV calling pipeline

"""

from os.path import split as pathsplit
from argparse import ArgumentParser

ARGPARSER = ArgumentParser(description='Argument for the SNV pipeline',
                                   prefix_chars='-')

ARGPARSER.add_argument('-project_name',
                       help='name of the project folder and where to find the fastq files (default: sample_test)',
                       default="sample_test",
                       metavar='str')

ARGPARSER.add_argument('-cell_type',
                       help=' (HUMAN/MOUSE) default: HUMAN',
                       default="HUMAN",
                       metavar='str')

ARGPARSER.add_argument('-read_length',
                       help=' star index read length (default: 51)',
                       default=51,
                       type=int,
                       metavar='int')

ARGPARSER.add_argument('-star_nb_threads',
                       help=' number of threads for STAR analysis (default 12)',
                       default=12,
                       type=int,
                       metavar='int')

ARGPARSER.add_argument('-snv_nb_threads',
                       help=' number of SNV calling pipelines executed in parallel (default 3)',
                       default=3,
                       type=int,
                       metavar='int')

ARGS = ARGPARSER.parse_known_args()[0]

# Project name. Used to create folder
PROJECT_NAME = ARGS.project_name
# type of the dataset (human or mouse). Used to select reference genomes
CELL_TYPE = ARGS.cell_type
# valid sequencing machine for picard tools:
# ILLUMINA, SLX, SOLEXA, SOLID, 454, LS454, COMPLETE, PACBIO,
# IONTORRE NT, CAPILLARY, HELICOS, UNKNOWN
PLATEFORM = 'ILLUMINA'
# Read length used to create star index for reference genome
STAR_INDEX_READ_LENGTH = ARGS.read_length

############ FOLDER ARCHITECTURE  ####################################
#Alias to define the GLOBAL_DATA_ROOT, OUTPUT_ROOT and PROG_ROOT
# (could be overloaded using reference paths)
USER = 'opoirion'
# Alias to define the root folder for reference data
# (could be overloaded using reference paths)
GLOBAL_DATA_ROOT = '/data/'
# Alias to define the output folder
OUTPUT_ROOT = '/data/results/'
# Alias to define the folder containing softwares.
# (could be overloaded using reference paths)
PROG_ROOT = '/prog/'
# Absolute path for the .soft file (dataset description) from NCBI
SOFT_PATH = "{0}/{1}/{1}.soft".format(GLOBAL_DATA_ROOT, PROJECT_NAME)
######################################################################

############ STANDART VARIABLE #######################################
TYPE_VAR = {
    'HUMAN': {
        # gtf file containing annotated human genes
        'ANNOTATION_PATH': "/data/Illumina_hg19/Annotation/genes.gtf",
        # folder which will contains the STAR index using human genome
        'STAR_INDEX_PATH': "{0}/Illumina_hg19/Sequences/STARindex".format(OUTPUT_ROOT),
        # folder which will contains the BSSEQ index using human genome
        'BSSEQ_INDEX_PATH': "/data/Illumina_hg19/Sequences/BSSEQindex".format(OUTPUT_ROOT),
        # human reference fasta (.fa) file
        'REF_GENOME': "/data/Illumina_hg19/Sequences/WholeGenomeFasta/genome.fa",
        # Reference human genome used
        'ORGANISM': 'hg19',
        # reference variant database used. The last version can be downloaded from:
        # ftp://ftp.ncbi.nih.gov/snp/organisms/ (human_9607_b{version}_p2)
        'DBSNP': "/data/Illumina_hg19/vcf/dbsnp_138.hg19.reduced.vcf",\
        'VCF_RESOURCES': [
            # Other reference variant resources.
            # Can be downloaded from ftp://ftp.broadinstitute.org/bundle/2.8/hg19
            # "/data/hg19/vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
            # Indel variant reference database
            # can be downloaded from ftp://ftp.broadinstitute.org/bundle/2.8/hg19
            # "/data/hg19/vcf/1000G_phase1.indels.hg19.sites.vcf",
            ]
    },
    'MOUSE': {
        # gtf file containing annotated mouse genes
        'ANNOTATION_PATH': "/data/Mus_musculus/UCSC/mm10/Annotation/genes.gtf",
        # folder which will contains the STAR index using mouse genome
        'STAR_INDEX_PATH': "{0}/Mus_musculus/UCSC/mm10/Sequence/STARindex".format(OUTPUT_ROOT),
        # folder which will contains the BS-SEQ index using mouse genome
        'BSSEQ_INDEX_PATH': "{0}/Mus_musculus/UCSC/mm10/Sequence/BSSEQindex".format(OUTPUT_ROOT),
        # Mouse reference fasta (.fa) file
        'REF_GENOME': "/data/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa",
        # Reference mouse genome used
        'ORGANISM': 'mm10',
        # reference variant database used. This version can be downloaded from:
        # ftp://ftp-mouse.sanger.ac.uk/REL-1303- SNPs_Indels-GRCm38/.
        'DBSNP': "/data/Mus_musculus/UCSC/mm10/vcf/mgp.v3.snps.rsIDdbSNPv137_ordered.reduced.vcf",
        # reference indel variant database used. This version can be downloaded from:
        # ftp://ftp-mouse.sanger.ac.uk/REL-1303- SNPs_Indels-GRCm38/.
        # Mouse VCF files must be sorted toward the sequence dictionnary of the mouse reference genome using SortVCF function from picard-tools
        'VCF_RESOURCES': [
            # "/data/mm10/vcf/mgp.v3.indels.rsIDdbSNPv137_ordered.vcf"
            ]
    }
}
######################################################################

############ MOUSE/ HUMAN ############################################
REF_GENOME = TYPE_VAR[CELL_TYPE]['REF_GENOME']
ANNOTATION_PATH = TYPE_VAR[CELL_TYPE]['ANNOTATION_PATH']
STAR_INDEX_PATH = TYPE_VAR[CELL_TYPE]['STAR_INDEX_PATH']
ORGANISM = TYPE_VAR[CELL_TYPE]['ORGANISM']
DBSNP = TYPE_VAR[CELL_TYPE]['DBSNP']
VCF_RESOURCES = TYPE_VAR[CELL_TYPE]['VCF_RESOURCES']
BSSEQ_INDEX_PATH = TYPE_VAR[CELL_TYPE]['BSSEQ_INDEX_PATH']
######################################################################

############# DATASET ################################################
# Absolute path for fastq files.
# Fastq files must be organised using one folder for one SRX experiment
FASTQ_PATH = "{0}/{1}/fastq/".format(OUTPUT_ROOT, PROJECT_NAME)
# output path
PATH_OUTPUT = "{0}/{1}/".format(OUTPUT_ROOT, PROJECT_NAME)
#specific string pattern that a folder name must match
SPECIFIC_FILENAME_PATTERN = ""
######################################################################

################ RRBS reads specific PREPROCESSING ###################
# Used aligner (star for reads from gene expression
# bismark / BS-seeker2 (RRBS read alignment))
USED_ALIGNER = 'STAR'
# Are the reads from the bislufite pipeline for SNV calling?
ARE_READS_BISULFITE = False
# specific trimming preprocessing for RRBS reads
DO_TRIMGALORE = True
######################################################################

######################## MONOVAR #####################################
MONOVAR_REP = '{0}/monovar/'.format(PROG_ROOT)
MONOVAR_SAMTOOLS = '{0}/external/samtools/samtools'.format(MONOVAR_REP)
######################################################################

############# SOFTWARE ###############################################
# Available java version. Must be > 1.8
JAVA = "java"
# Max memory used by Java
JAVA_MEM = "-Xmx110g"
# GATK folder where can be found GATK software
GATK_DIR = "{0}/GATK/".format(PROG_ROOT)
# GATK jar name
GATK_JAR = "GenomeAnalysisTK.jar"
# picard-tools software
PICARD_DIR = "{0}/picard-tools-2.1.1/".format(PROG_ROOT)
# Perl
PERL = 'perl'
# python
PYTHON = 'python'
#BOWTIE ALIGNER (for BSSEEKER and bismark)
BOWTIE_REP = '/usr/bin/'
# software for RRBS bisulfite reads preprocessing
TRIMGALORE_REP = '{0}/TrimGalore/'.format(PROG_ROOT)
# BSseeker2 software to call methylation reads
BSSEEKER2_REP = '{0}/BSseeker2/'.format(PROG_ROOT)
# BS-Snper (SNP calling for bisulfite reads)
BSSNPER = '{0}/BS-Snper/BS-Snper.pl'.format(PROG_ROOT)
# bismark software for RRBS alignment
BISMARK_SOFTWARE = '{0}/Bismark/bismark'.format(PROG_ROOT)
# STAR aligner software
PATH_STAR_SOFTWARE = "{0}/STAR/bin/Linux_x86_64_static/STAR"\
                          .format(PROG_ROOT)
# fastqc software [OPTIONAL]
FASTQC = "fastqc"
# snpEff software (vcf annotation) [OPTIONAL]
SNPEFF = '{0}/snpEff/snpEff.jar'.format(PROG_ROOT)
# required snpEff databases (vcf annotation) [OPTIONAL]
SNPEFF_DICT = {'MOUSE': 'GRCm38.82',
               'HUMAN': 'GRCh37.75'}
SNPEFF_DB = SNPEFF_DICT[CELL_TYPE]
# SAMtools
SAMTOOLS = '{0}/samtools-1.5/bin/samtools'.format(PROG_ROOT)
######################################################################

#############  STAR #################################################
# Number of threads used when using STAR aligner
STAR_THREADS = ARGS.star_nb_threads
# output path for STAR results
OUTPUT_PATH_STAR = PATH_OUTPUT + "/star/"
#####################################################################

############ SNV CALLING PIPELINE ###################################
# output path for SNVs inferred
OUTPUT_PATH_GATK =  PATH_OUTPUT + '/snv_pipeline_GATK/'
# Number of SNV calling processes launched in parallel
NB_PROCESS_SNV = ARGS.snv_nb_threads
####################################################################

########### FREEBAYES SNV CALLING PIPELINE ##########################
OUTPUT_PATH_FREEBAYES = PATH_OUTPUT + '/snv_pipeline_freebayes/'
PATH_OPOSSUM = '{0}/Opossum/'.format(PROG_ROOT)
PATH_FREEBAYES = '{0}/freebayes/bin/freebayes'.format(PROG_ROOT)
####################################################################

############ COMPUTE DISTANCE MATRIX [OPTIONAL] ##############################
# software to infer gene expressions count with raw count for each single cells
# [OPTIONAL]
FEATURE_COUNT = "featureCounts"
# path for gene expression matrices [OPTIONAL]
MATRIX_OUTPUT_PATH = "{0}/{1}/expression_profile/"\
                     .format(OUTPUT_ROOT, PROJECT_NAME)
###############################################################################


######################## SNV SIMULATION #######################################
SIMULATED_REF_GENOME = None

if SIMULATED_REF_GENOME:
    MUTATION_FILE = '{0}/Simulated{1}Mut/sim_snv.bed'.format(
        pathsplit(pathsplit(REF_GENOME)[0])[0], SIMULATED_REF_GENOME)
    SEQUENCES_PATH = pathsplit(pathsplit(REF_GENOME)[0])[0]
    SIM_GENOME_DIR = '{0}/Simulated{1}Mut/'.format(SEQUENCES_PATH, SIMULATED_REF_GENOME)
    REF_GENOME_ORIGINAL = REF_GENOME[:]
    REF_GENOME = '{0}/genome.fa'.format(SIM_GENOME_DIR)
else:
    MUTATION_FILE = None
    SEQUENCES_PATH = None
    SIM_GENOME_DIR = None
    REF_GENOME_ORIGINAL = None
###############################################################################
