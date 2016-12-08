#! /usr/bin/python

""" SNV calling pipeline"""

from os import popen
from os.path import isdir
from os.path import isfile
from os.path import getsize
from subprocess import Popen
from subprocess import PIPE

from distutils.dir_util import mkpath
from shutil import rmtree
from shutil import copyfile
from shutil import move
from sys import stdout as STDOUT
from sys import argv
from random import randint
from random import random
from time import sleep
from time import time

from garmire_SNV_calling.config import OUTPUT_PATH_SNV as CURRENT_PATH
from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import JAVA_MEM
from garmire_SNV_calling.config import PICARD_DIR
from garmire_SNV_calling.config import GATK_DIR
from garmire_SNV_calling.config import PLATEFORM
from garmire_SNV_calling.config import ORGANISM
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import DBSNP
from garmire_SNV_calling.config import VCF_RESOURCES


############ VARIABLES ################
SRR_TO_PROCESS = "" # for debug purpose
#######################################


class ProcessSNVCalling():
    """ """
    def __init__(self,
                 current_path=CURRENT_PATH,
                 path_to_data=PATH_OUTPUT,
                 id="1",
                 clean_tmp=True,
    ):
        self.current_path = current_path
        self.path_to_data = path_to_data
        self.id = str(id)
        self.stdout = None
        self.tmppath = None
        self.time_start = None
        self.clean_tmp = clean_tmp
        self.srr_to_process = None

    def process(self, srr_to_process=SRR_TO_PROCESS):
        """
        process one star bam file with snv calling pipeline
        """
        self._init_process(srr_to_process)
        self._launch_picard_readgroups()
        self._launch_picard_markduplicates()
        self._launch_gatk_cigar()
        self._launch_gatk_realigner_target_creator()
        self._launch_gatk_realigner_indel()
        self._launch_gatk_base_recalibrator()
        self._launch_gatk_print_reads()
        self._launch_gatk_variant_calling()
        self._launch_gatk_variant_filtering()
        self._finish_process()

    def _init_process(self, srr_to_process):
        """mk tmp folders... """
        self.time_start = time()
        self.srr_to_process = srr_to_process
        self.tmppath = self.current_path + "/tmp/" + self.id
        self.respath = self.current_path + "/data/" + srr_to_process

        if isfile(self.respath + '/snv_filtered.vcf') \
                and getsize(self.respath + '/snv_filtered.vcf'):
            err = 'error file : {0} already exists!'\
                .format(self.respath + '/snv_filtered.vcf')
            raise Exception(err)

        sleep(2 * random())
        if not isdir(self.tmppath):
            mkpath(self.tmppath)

        if self.clean_tmp:
            popen("rm {0}/*".format(self.tmppath)).read()

        self.stdout = open(self.tmppath + '/stdout.log', 'a+')

        self.stdout.write('\n\n######## file id {0} ########\n'\
                          .format(srr_to_process))
        star_path = self.path_to_data + "/star/"  + srr_to_process

        if not isfile(star_path + '/Aligned.sortedByCoord.out.bam')\
                or not getsize(star_path + '/Aligned.sortedByCoord.out.bam'):
            err = 'error file : {0} not found or empty!'\
                .format(star_path + '/Aligned.sortedByCoord.out.bam')
            raise Exception(err)

        copyfile("{0}/Aligned.sortedByCoord.out.bam".format(star_path),
                 "{0}/Aligned.sortedByCoord.out.bam".format(self.tmppath))


    def _finish_process(self):
        """mk res folders... """

        if not isdir(self.respath):
            mkpath(self.respath)

        copyfile(self.tmppath + '/snv_filtered.vcf',
                 self.respath + '/snv_filtered.vcf')
        copyfile(self.tmppath + '/snv_filtered.vcf.idx',
                 self.respath  + '/snv_filtered.vcf.idx')
        copyfile(self.tmppath + '/stdout.log',
                 self.respath + '/stdout.log')
        self._run_cmd('echo "#### FINISHED ####'\
                      ' ALL PROCESS DONE FOR: {0} in {1} s"'\
                      .format(self.srr_to_process, time() - self.time_start))

    def _launch_picard_readgroups(self):
        """
        launch picard AddOrReplaceReadGroups
        """
        popen("rm {0}/rg_added_sorted.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING PICARD READGROUPS ########\n"')

        cmd = "{0} {1} -jar {2}/picard.jar AddOrReplaceReadGroups" \
              " I={3}/Aligned.sortedByCoord.out.bam"\
              " O={3}/rg_added_sorted.bam" \
              " SO=coordinate" \
              " RGID={4}" \
              " RGPU={4}" \
              " RGSM={4}" \
              " RGPL={5}" \
              " RGLB={6}" \
              .format(JAVA,
                      JAVA_MEM,
                      PICARD_DIR,
                      self.tmppath,
                      self.id,
                      PLATEFORM,
                      ORGANISM
              )
        self._run_cmd(cmd)

    def _launch_picard_markduplicates(self):
        """
        launch picard MarkDuplicates
        """
        popen("rm {0}/dedupped.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING PICARD MARKDUPLICATES ########\n"')

        cmd = "{0} {1} -jar {2}/picard.jar MarkDuplicates" \
              " I={3}/rg_added_sorted.bam"\
              " O={3}/dedupped.bam" \
              " M={3}/output.metrics" \
              " CREATE_INDEX=true" \
              " VALIDATION_STRINGENCY=SILENT" \
              .format(JAVA,
                      JAVA_MEM,
                      PICARD_DIR,
                      self.tmppath,
              )
        self._run_cmd(cmd)

    def _launch_picard_reordersam(self):
        """
        launch picard REORDERSAM
        """
        popen("rm {0}/reordered.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING PICARD REORDERSAM ########\n"')

        cmd = "{0} {1} -jar {2}/picard.jar ReorderSam" \
              " I={3}/dedupped.bam" \
              " O={3}/dedupped.bam" \
              " R={4}"\
              " CREATE_INDEX=TRUE" \
              .format(JAVA,
                      JAVA_MEM,
                      PICARD_DIR,
                      self.tmppath,
                      REF_GENOME
              )
        self._run_cmd(cmd)

    def _launch_gatk_cigar(self):
        """
        Running cigar string split and mapq 255 fix GATK
        """
        popen("rm {0}/split.ba*".format(self.tmppath)).read()
        self._run_cmd('echo "\n\n######## LAUNCHING CIGAR ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T SplitNCigarReads" \
        " -I {3}/dedupped.bam" \
        " -o {3}/split.bam" \
        " -R {4}" \
        " -rf ReassignOneMappingQuality" \
        " -RMQF 255" \
        " -RMQT 60" \
        " -U ALLOW_N_CIGAR_READS" \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME
        )

        self._run_cmd_fix_quality(cmd, to_rm='split.ba*')

    def _launch_gatk_realigner_target_creator(self):
        """
        Running Realignment Target creator
        """
        popen("rm {0}/forRealigner.intervals".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING REALIGNER TARGET CREATOR ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T RealignerTargetCreator" \
        " -I {3}/split.bam" \
        " -o {3}/forRealigner.intervals"\
        " -R {4}" \
        " -nt 20 " \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME
              )

        for vcf in VCF_RESOURCES:
            cmd += " -known {0}".format(vcf)

        self._run_cmd_fix_quality(cmd, to_rm='forRealigner.intervals', resolve='hard')

    def _launch_gatk_realigner_indel(self):
        """
        Running Realignment
        """
        popen("rm {0}/realigned.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING REALIGNER INDEL ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T IndelRealigner" \
        " -I {3}/split.bam" \
        " -targetIntervals {3}/forRealigner.intervals"\
        " --out {3}/realigned.bam" \
        " -R {4}" \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME,
              )

        for vcf in VCF_RESOURCES:
            cmd += " -known {0}".format(vcf)

        self._run_cmd(cmd)

    def _launch_gatk_base_recalibrator(self):
        """
        Running base recalibration
        """
        popen("rm {0}/recal_data.csv".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING RECALIBRATION STEP 1 ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T BaseRecalibrator" \
        " -I {3}/realigned.bam" \
        " -o {3}/recal_data.csv" \
        " -R {4}" \
        " -nct 20" \
        " --knownSites {5}" \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME,
                DBSNP
              )

        for vcf in VCF_RESOURCES:
            cmd += " --knownSites {0}".format(vcf)

        self._run_cmd_fix_quality(cmd, to_rm='recal_data.csv', resolve='hard')

    def _launch_gatk_print_reads(self):
        """
        Running base recalibration STEP 2
        """
        popen("rm {0}/recal.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING RECALIBRATION STEP 2 ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T PrintReads" \
        " -I {3}/realigned.bam" \
        " --out {3}/recal.bam" \
        " -R {4}" \
        " -BQSR {3}/recal_data.csv" \
        " -nct 20" \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME
              )

        self._run_cmd_fix_quality(cmd, to_rm='recal.bam', resolve='hard')

    def _launch_gatk_variant_calling(self):
        """
        variant calling
        """
        popen("rm {0}/snv_raw.vcf".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING VARIANT CALLING ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T HaplotypeCaller" \
        " -I {3}/recal.bam" \
        " -o {3}/snv_raw.vcf" \
        " -R {4}" \
        " --dbsnp {5}" \
        " -dontUseSoftClippedBases" \
        " -stand_call_conf 20.0" \
        " -stand_emit_conf 20.0" \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME,
                DBSNP
              )

        self._run_cmd(cmd)

    def _launch_gatk_variant_filtering(self):
        """
        variant filtering
        """
        popen("rm {0}/snv_filtered.vcf".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n######## LAUNCHING VARIANT FILTERING ########\n"')

        cmd = "{0} {1} -jar {2}/GenomeAnalysisTK.jar -T VariantFiltration" \
        " -V {3}/snv_raw.vcf" \
        " -o {3}/snv_filtered.vcf" \
        " -R {4}" \
        " -cluster 3" \
        " -filterName FS" \
        ' -filter "FS > 30.0"' \
        " -filterName QD" \
        ' -filter "QD < 2.0"' \
        .format(JAVA,
                JAVA_MEM,
                GATK_DIR,
                self.tmppath,
                REF_GENOME,
                DBSNP
              )

        self._run_cmd(cmd)

    def _run_cmd(self, cmd):
        """run cmd"""
        stdout_read = open(self.tmppath + '/stdout.log', 'r')
        stdout_read.seek(0, 2)

        process = Popen(cmd,
                        stdout=PIPE,
                        stderr=PIPE,
                        shell=True)

        c = process.stdout.read(1)
        e = process.stderr.read(1)

        while process.poll() == None or c or e:

            STDOUT.write(c)
            self.stdout.write(c)
            STDOUT.write(e)
            self.stdout.write(e)
            STDOUT.flush()
            self.stdout.flush()

            c = process.stdout.read(1)
            e = process.stderr.read(1)

        process.communicate()

        if process.returncode:
            raise Exception('{0} raise non 0 return code!\n'\
                            .format(cmd))

    def _run_cmd_fix_quality(self, cmd, to_rm, resolve='soft'):
        """ """
        try:
            self._run_cmd(cmd)
        except Exception as e:
            self._run_cmd('echo "\n\nERROR DETECTED.' \
                          'Try correcting missencoded quality score"')

            if resolve == 'hard':
                cmd += ' --allow_potentially_misencoded_quality_scores'
            else:
                cmd += ' --fix_misencoded_quality_scores'

            popen("rm {0}/{1}".format(self.tmppath, to_rm)).read()
            self._run_cmd(cmd)
