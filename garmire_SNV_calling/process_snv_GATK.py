#! /usr/bin/python

""" process one SSR with GATK pipeline SNV"""

from os import popen
from os.path import isdir
from os.path import isfile
from os.path import getsize
from subprocess import Popen
from subprocess import PIPE

from distutils.dir_util import mkpath

from shutil import copyfile

from sys import stdout as STDOUT
from sys import argv
from random import randint
from random import random
from time import sleep
from time import time

from glob import glob

from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import JAVA_MEM
from garmire_SNV_calling.config import PICARD_DIR
from garmire_SNV_calling.config import GATK_DIR
from garmire_SNV_calling.config import GATK_JAR


############ VARIABLES ############################################
SRR_TO_PROCESS = "" # for debug purpose
PROCESS_ID = randint(0, 1000000)

if "--specific_folder" in argv:
    SRR_TO_PROCESS = argv[
        argv.index("--specific_folder") + 1]
if "--process_id" in argv:
    PROCESS_ID = int(argv[
        argv.index("--process_id") + 1])

if "--path_to_data" in argv:
    PATH_TO_DATA = argv[
        argv.index("--path_to_data") + 1]
    OUTPUT_PATH =  PATH_TO_DATA + '/snv_pipeline_raw/'
else:
    from garmire_SNV_calling.config import PATH_TO_DATA
    from garmire_SNV_calling.config import OUTPUT_PATH_GATK as OUTPUT_PATH

if "--plateform" in argv:
    PLATEFORM = argv[
        argv.index("--plateform") + 1]
else:
    from garmire_SNV_calling.config import PLATEFORM

if "--organism" in argv:
    ORGANISM = argv[
        argv.index("--organism") + 1]
else:
    from garmire_SNV_calling.config import ORGANISM

if "--ref_genome" in argv:
    REF_GENOME = argv[
        argv.index("--ref_genome") + 1]
else:
    from garmire_SNV_calling.config import REF_GENOME

if "--dbsnp" in argv:
    DBSNP = argv[
        argv.index("--dbsnp") + 1]
else:
    from garmire_SNV_calling.config import DBSNP

if "--vcf_resources" in argv:
    VCF_RESOURCES = eval(argv[
        argv.index("--vcf_resources") + 1])
else:
    from garmire_SNV_calling.config import VCF_RESOURCES


if "--ignore_already_exists" in argv:
    IGNORE_ALREADY_EXISTS = True
else:
    IGNORE_ALREADY_EXISTS = False

if "--clean_tmp" in argv:
    CLEAN_TMP = True
else:
    CLEAN_TMP = False


###################################################################


def main():
    process_GATK_snv = ProcessGATKSNV(id=PROCESS_ID)
    process_GATK_snv.process()


class ProcessGATKSNV():
    """ """
    def __init__(self,
                 bam_file_name='',
                 srr_to_process=SRR_TO_PROCESS,
                 output_path=OUTPUT_PATH,
                 path_to_data=PATH_TO_DATA,
                 java=JAVA,
                 java_mem=JAVA_MEM,
                 picard_dir=PICARD_DIR,
                 gatk_dir=GATK_DIR,
                 plateform=PLATEFORM,
                 organism=ORGANISM,
                 ref_genome=REF_GENOME,
                 dbsnp=DBSNP,
                 vcf_resources=VCF_RESOURCES,
                 gatk_jar=GATK_JAR,
                 id="1",
                 ignore_already_exists=IGNORE_ALREADY_EXISTS,
                 clean_tmp=CLEAN_TMP,
                 respath=None,
    ):

        self.respath = respath

        self.output_path = output_path
        self.path_to_data = path_to_data
        self.srr_to_process = srr_to_process

        self.bam_file_name = bam_file_name

        self.java = java
        self.java_mem = java_mem
        self.picard_dir = picard_dir
        self.gatk_dir = gatk_dir
        self.plateform = plateform
        self.organism = organism[:]
        self.ignore_already_exists = ignore_already_exists

        if self.organism == 'HUMAN':
            self.organism = 'hg19'

        elif self.organism == 'MOUSE':
            self.organism = 'mm10'

        self.ref_genome = ref_genome
        self.dbsnp = dbsnp
        self.vcf_resources = vcf_resources

        self.id = str(id)
        self.stdout = None
        self.tmppath = None
        self.time_start = None
        self.bam_file_path = None
        self.clean_tmp = clean_tmp

    def process(self, srr_to_process=None):
        """
        process one star bam file with snv calling pipeline
        """
        if srr_to_process:
            self.srr_to_process = srr_to_process

        msg = self._init_process()

        if msg:
            print(msg)
            self.stdout.write(msg)
            return

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
        self._rm_tmp_file()

    def process_exome(self, srr_to_process=None):
        """
        process one star bam file with snv calling pipeline
        """
        if srr_to_process:
            self.srr_to_process = srr_to_process


        msg = self._init_process()

        if msg:
            print(msg)
            self.stdout.write(msg)
            return

        self._launch_picard_readgroups()
        self._launch_picard_buildbamindex(name='rg_added_sorted')
        self._launch_picard_markduplicates()
        self._launch_gatk_base_recalibrator(input_name='dedupped')
        self._launch_gatk_print_reads(input_name='dedupped')
        self._launch_gatk_variant_calling()
        self._launch_gatk_variant_filtering()
        self._finish_process()
        self._rm_tmp_file()

    def _init_process(self):
        """mk tmp folders... """
        self.time_start = time()
        self.tmppath = self.output_path + "/tmp/" + self.id

        if not self.respath:
            self.respath = self.output_path + \
                           "/data/" + self.srr_to_process

        sleep(2 * random())
        if not isdir(self.tmppath):
            mkpath(self.tmppath)

        if self.clean_tmp and glob('{0}/*'.format(self.tmppath)):
            popen("rm {0}/*".format(self.tmppath)).read()

        self.stdout = open(self.tmppath + '/stdout.log', 'a+')

        self.stdout.write('\n\n######## file id {0} ########\n'\
                          .format(self.srr_to_process))

        if isfile(self.respath + '/snv_filtered.vcf') \
                and getsize(self.respath + '/snv_filtered.vcf'):
            msg = 'file : {0} already exists!'\
                .format(self.respath + '/snv_filtered.vcf')
            print(msg)

            if self.ignore_already_exists:
                print('continuing anyway...')
            else:
                return msg

        if not self.bam_file_name:
            self.bam_file_path = '{0}{1}/Aligned.sortedByCoord.out.bam'.format(
                self.path_to_data +  "/star/" ,  self.srr_to_process)
        else:
            self.bam_file_path = self.path_to_data + self.bam_file_name

        if not isfile(self.bam_file_path)\
                or not getsize(self.bam_file_path):
            err = 'error file : {0} not found or empty!'\
                .format(self.bam_file_path)
            raise Exception(err)

        copyfile("{0}".format(self.bam_file_path),
                 "{0}/Aligned.sortedByCoord.out.bam".format(self.tmppath))

    def _finish_process(self):
        """mk res folders... """

        if not isdir(self.respath):
            mkpath(self.respath)

        copyfile(self.tmppath + '/snv_filtered.vcf',
                 self.respath + '/snv_filtered_freebayes.vcf')

        if isfile(self.tmppath + '/snv_filtered.vcf.idx'):
            copyfile(self.tmppath + '/snv_filtered.vcf.idx',
                     self.respath  + '/snv_filtered_freebayes.vcf.idx')

        self.stdout.write('''\n #### FINISHED #### \n
ALL PROCESS DONE FOR: {0} in {1} s
        '''.format(self.srr_to_process, time() - self.time_start))

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
              .format(self.java,
                      self.java_mem,
                      self.picard_dir,
                      self.tmppath,
                      self.id,
                      self.plateform,
                      self.organism
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
              .format(self.java,
                      self.java_mem,
                      self.picard_dir,
                      self.tmppath,
              )
        self._run_cmd(cmd)

    def _launch_picard_buildbamindex(self, name='dedupped'):
        """
        launch picard buildbamindex
        """
        popen("rm {0}/{1}.bai".format(self.tmppath, name)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING PICARD BuildBamIndex ########\n"')

        cmd = "{0} {1} -jar {2}/picard.jar BuildBamIndex" \
              " I={3}/{4}.bam" \
              " TMP_DIR={3}" \
              .format(self.java,
                      self.java_mem,
                      self.picard_dir,
                      self.tmppath,
                      name,
              )
        self._run_cmd(cmd)

    def _launch_picard_sortsam(self):
        """
        launch picard SORTSAM
        """
        popen("rm {0}/sorted.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING PICARD REORDERSAM ########\n"')

        cmd = "{0} {1} -jar {2}/picard.jar SortSam" \
              " I={3}/dedupped.bam" \
              " O={3}/sorted.bam" \
              " SORT_ORDER=coordinate" \
              " TMP_DIR={3}" \
              " CREATE_INDEX=TRUE" \
              .format(self.java,
                      self.java_mem,
                      self.picard_dir,
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
              " O={3}/dedupped_reodered.bam" \
              " R={4}"\
              " CREATE_INDEX=TRUE" \
              .format(self.java,
                      self.java_mem,
                      self.picard_dir,
                      self.tmppath,
                      self.ref_genome
              )
        self._run_cmd(cmd)

    def _launch_gatk_cigar(self):
        """
        Running cigar string split and mapq 255 fix GATK
        """
        popen("rm {0}/split.ba*".format(self.tmppath)).read()
        self._run_cmd('echo "\n\n######## LAUNCHING CIGAR ########\n"')

        cmd = "{0} {1} -jar {2}/{5} -T SplitNCigarReads" \
        " -I {3}/dedupped.bam" \
        " -o {3}/split.bam" \
        " -R {4}" \
        " -rf ReassignOneMappingQuality" \
        " -RMQF 255" \
        " -RMQT 60" \
        " -U ALLOW_N_CIGAR_READS" \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                self.gatk_jar
        )

        self._run_cmd_fix_quality(cmd, to_rm='split.ba*')

    def _launch_gatk_realigner_target_creator(self, input_name='split.bam', resolve='hard'):
        """
        Running Realignment Target creator
        """
        popen("rm {0}/forRealigner.intervals".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING REALIGNER TARGET CREATOR ########\n"')

        cmd = "{0} {1} -jar {2}/{6} -T RealignerTargetCreator" \
        " -I {3}/{5}" \
        " -o {3}/forRealigner.intervals"\
        " -R {4}" \
        " -nt 20 " \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                input_name,
                self.gatk_jar
              )

        for vcf in self.vcf_resources:
            cmd += " -known {0}".format(vcf)

        self._run_cmd_fix_quality(cmd, to_rm='forRealigner.intervals', resolve=resolve)

    def _launch_gatk_realigner_indel(self):
        """
        Running Realignment
        """
        popen("rm {0}/realigned.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING REALIGNER INDEL ########\n"')

        cmd = "{0} {1} -jar {2}/{5} -T IndelRealigner" \
        " -I {3}/split.bam" \
        " -targetIntervals {3}/forRealigner.intervals"\
        " --out {3}/realigned.bam" \
        " -R {4}" \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                self.gatk_jar
              )

        for vcf in self.vcf_resources:
            cmd += " -known {0}".format(vcf)

        self._run_cmd(cmd)

    def _launch_gatk_base_recalibrator(self, input_name='realigned'):
        """
        Running base recalibration
        """
        popen("rm {0}/recal_data.csv".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING RECALIBRATION STEP 1 ########\n"')

        cmd = "{0} {1} -jar {2}/{6} -T BaseRecalibrator" \
        " -I {3}/{6}.bam" \
        " -o {3}/recal_data.csv" \
        " -R {4}" \
        " -nct 20" \
        " --knownSites {5}" \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                self.dbsnp,
                input_name,
                self.gatk_jar
              )

        for vcf in self.vcf_resources:
            cmd += " --knownSites {0}".format(vcf)

        self._run_cmd_fix_quality(cmd, to_rm='recal_data.csv', resolve='hard')

    def _launch_gatk_print_reads(self, input_name='realigned'):
        """
        Running base recalibration STEP 2
        """
        popen("rm {0}/recal.bam".format(self.tmppath)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING RECALIBRATION STEP 2 ########\n"')

        cmd = "{0} {1} -jar {2}/{6} -T PrintReads" \
        " -I {3}/{5}.bam" \
        " --out {3}/recal.bam" \
        " -R {4}" \
        " -BQSR {3}/recal_data.csv" \
        " -nct 20" \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                input_name,
                self.gatk_jar
              )

        self._run_cmd_fix_quality(cmd, to_rm='recal.bam', resolve='hard')

    def _launch_gatk_variant_calling(self, output_name='snv_raw.vcf'):
        """
        variant calling
        """
        popen("rm {0}/{1}".format(self.tmppath, output_name)).read()
        self._run_cmd(
            'echo "\n\n######## LAUNCHING VARIANT CALLING ########\n"')

        start_time = time()

        cmd = "{0} {1} -jar {2}/{7} -T HaplotypeCaller" \
        " -I {3}/recal.bam" \
        " -o {3}/{6}" \
        " -R {4}" \
        " --dbsnp {5}" \
        " -dontUseSoftClippedBases" \
        " -stand_call_conf 20.0" \
        " -stand_emit_conf 20.0" \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                self.dbsnp,
                output_name,
                self.gatk_jar

              )

        self._run_cmd(cmd)

        self._run_cmd(
            'echo "\n## GATK variant calling done in {0} s##\n"'.format(
                time() - start_time))

    def _launch_gatk_variant_filtering(
            self,
            input_name='snv_raw.vcf',
            output_name='snv_filtered.vcf'):
        """
        variant filtering
        """
        popen("rm {0}/{1}".format(self.tmppath, output_name)).read()
        self._run_cmd(
            'echo "\n######## LAUNCHING VARIANT FILTERING ########\n"')

        start_time = time()

        cmd = "{0} {1} -jar {2}/{7} -T VariantFiltration" \
        " -V {3}/{5}" \
        " -o {3}/{6}" \
        " -R {4}" \
        " -cluster 3" \
        " -filterName FS" \
        ' -filter "FS > 30.0"' \
        " -filterName QD" \
        ' -filter "QD < 2.0"' \
        .format(self.java,
                self.java_mem,
                self.gatk_dir,
                self.tmppath,
                self.ref_genome,
                input_name,
                output_name,
                self.gatk_jar
              )

        self._run_cmd(cmd)

        self._run_cmd(
            'echo "\n## GATK variant filtering done in {0} s##\n"'.format(
                time() - start_time))

    def _rm_tmp_file(self):
        """
        """
        if isdir(self.tmppath) and self.clean_tmp:
            cmd = 'rm -rf {0}'.format(self.tmppath)
            try:
                self._run_cmd(cmd)
            except Exception as e:
                print('#### error while trying to remove the tmp folder: {0}'\
                  .format(e))

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
        except Exception:
            self._run_cmd('echo "\n\nERROR DETECTED.' \
                          'Try correcting missencoded quality score"')

            if resolve == 'hard':
                cmd += ' --allow_potentially_misencoded_quality_scores'
            else:
                cmd += ' --fix_misencoded_quality_scores'

            popen("rm {0}/{1}".format(self.tmppath, to_rm)).read()
            self._run_cmd(cmd)

if __name__ == "__main__":
    main()
