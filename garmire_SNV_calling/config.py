"""
config file for SNV calling pipeline

"""

PROJECT_NAME = 'jones_pancreatic_cancer'
# Project name. Used to create folder
CELL_TYPE = 'HUMAN'
# type of the dataset (human or mouse). Used to select reference genomes
PLATEFORM = 'ILLUMINA'
# valid sequencing machine for picard tools:
# ILLUMINA, SLX, SOLEXA, SOLID, 454, LS454, COMPLETE, PACBIO,
# IONTORRE NT, CAPILLARY, HELICOS, UNKNOWN
STAR_INDEX_READ_LENGTH = 50
# Read length used to create star index for reference genome

############ FOLDER ARCHITECTURE  ####################################
USER = 'opoirion'
#Alias to define the GLOBAL_DATA_ROOT, OUTPUT_ROOT and PROG_ROOT
# (could be overloaded using reference paths)
GLOBAL_DATA_ROOT = '/data/{0}/'.format(USER)
# Alias to define the root folder for reference data
# (could be overloaded using reference paths)
OUTPUT_ROOT = '/home/{0}/data/'.format(USER)
# Alias to define the output folder
PROG_ROOT = '/home/{0}/prog/'.format(USER)
# Alias to define the folder containing softwares.
# (could be overloaded using reference paths)
SOFT_PATH = "{0}/{1}/{1}.soft".format(GLOBAL_DATA_ROOT, PROJECT_NAME)
# Absolute path for the .soft file (dataset description) from NCBI
######################################################################

############ STANDART VARIABLE #######################################
TYPE_VAR = {
    'HUMAN': {
        'ANNOTATION_PATH': "{0}/Illumina_hg19/Annotation/genes.gtf"\
        .format(GLOBAL_DATA_ROOT),
        # gtf file containing annotated human genes
        'STAR_INDEX_PATH': "{0}/Illumina_hg19/Sequences/STARindex"\
        .format(GLOBAL_DATA_ROOT),
        # folder which will contains the STAR index using human genome
        'REF_GENOME': "{0}/Illumina_hg19/Sequences/WholeGenomeFasta/genome.fa"\
        .format(GLOBAL_DATA_ROOT),
        # human reference fasta (.fa) file
        'ORGANISM': 'hg19',
        # Reference human genome used
        'DBSNP': "{0}/Illumina_hg19/vcf/dbsnp_138.hg19.vcf"\
        # reference variant database used. The last version can be downloaded from:
        # ftp://ftp.ncbi.nih.gov/snp/organisms/ (human_9607_b{version}_p2)
        .format(GLOBAL_DATA_ROOT),
        'VCF_RESOURCES': [
            "{0}/Illumina_hg19/vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"\
            .format(GLOBAL_DATA_ROOT),
            # Other reference variant resources.
            # Can be downloaded from ftp://ftp.broadinstitute.org/bundle/2.8/hg19
            "{0}/Illumina_hg19/vcf/1000G_phase1.indels.hg19.sites.vcf"\
            .format(GLOBAL_DATA_ROOT),
            # Indel variant reference database
            # can be downloaded from ftp://ftp.broadinstitute.org/bundle/2.8/hg19
            ]
    },
    'MOUSE': {
        'ANNOTATION_PATH': "{0}/Mus_musculus/UCSC/mm10/Annotation/genes.gtf"\
        .format(GLOBAL_DATA_ROOT),
        # gtf file containing annotated mouse genes
        'STAR_INDEX_PATH': "{0}/Mus_musculus/UCSC/mm10/Sequence/STARindex"\
        .format(GLOBAL_DATA_ROOT),
        # folder which will contains the STAR index using mouse genome
        'REF_GENOME': "{0}/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"\
        .format(GLOBAL_DATA_ROOT),
        # Mouse reference fasta (.fa) file
        'ORGANISM': 'mm10',
        # Reference mouse genome used
        'DBSNP': "{0}/Mus_musculus/vcf/mgp.v3.snps.rsIDdbSNPv137_ordered.vcf"\
        .format(GLOBAL_DATA_ROOT),
        # reference variant database used. This version can be downloaded from:
        # ftp://ftp-mouse.sanger.ac.uk/REL-1303- SNPs_Indels-GRCm38/.
        'VCF_RESOURCES': [
            "{0}/Mus_musculus/vcf/mgp.v3.indels.rsIDdbSNPv137_ordered.vcf"\
            .format(GLOBAL_DATA_ROOT)
        # reference indel variant database used. This version can be downloaded from:
        # ftp://ftp-mouse.sanger.ac.uk/REL-1303- SNPs_Indels-GRCm38/.
        # Mouse VCF files must be sorted toward the sequence dictionnary of the mouse reference genome using SortVCF function from picard-tools
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
######################################################################

############# DATASET ###############################################
FASTQ_PATH = "{0}/{1}/fastq/".format(GLOBAL_DATA_ROOT, PROJECT_NAME)
# Absolute path for fastq files.
# Fastq files must be organised using one folder for one SRX experiment
PATH_OUTPUT = "{0}/{1}/".format(OUTPUT_ROOT, PROJECT_NAME)
# output path
SPECIFIC_FILENAME_PATTERN = ""
#####################################################################

############# SOFTWARE ###############################################
JAVA = "{0}/jdk1.8.0_77/bin/java".format(PROG_ROOT)
# Available java version. Must be > 1.8
JAVA_MEM = "-Xmx110g"
# Max memory used by Java
GATK_DIR = "{0}/GATK/".format(PROG_ROOT)
# GATK folder where can be found GATK software
PICARD_DIR = "{0}/picard-tools-2.1.1/".format(PROG_ROOT)
# picard-tools software
PATH_STAR_SOFTWARE = "{0}/STAR/bin/Linux_x86_64_static/STAR"\
                          .format(PROG_ROOT)
# STAR aligner software
FASTQC = "fastqc"
# fastqc software [OPTIONAL]
SNPEFF = '{0}/snpEff/snpEff.jar'.format(PROG_ROOT)
# snpEff software (vcf annotation) [OPTIONAL]
SNPEFF_DICT = {'MOUSE': 'GRCm38.82',
               'HUMAN': 'GRCh37.75'}
# required snpEff databases (vcf annotation) [OPTIONAL]
SNPEFF_DB = SNPEFF_DICT[CELL_TYPE]
######################################################################

#############  STAR #################################################
STAR_THREADS = 12
# Number of threads used when using STAR aligner
OUTPUT_PATH_STAR = PATH_OUTPUT + "/star/"
# output path for STAR results
#####################################################################

############ SNV CALLING PIPELINE ###################################
OUTPUT_PATH_SNV =  PATH_OUTPUT + '/snv_pipeline_results/'
# output path for SNVs inferred
NB_PROCESS_SNV = 2
# Number of SNV calling processes launched in parallel
####################################################################

############ COMPUTE DISTANCE MATRIX [OPTIONAL] ##############################
FEATURE_COUNT = "featureCounts"
# software to infer gene expressions count with raw count for each single cells
# [OPTIONAL]
MATRIX_OUTPUT_PATH = "{0}/{1}/expression_profile/"\
                     .format(OUTPUT_ROOT, PROJECT_NAME)
# path for gene expression matrices [OPTIONAL]
###############################################################################
