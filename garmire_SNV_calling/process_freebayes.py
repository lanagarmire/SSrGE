from config import SAMTOOLS
from config import PYTHON

from sys import argv
from os.path import isfile

from time import time

from garmire_SNV_calling.process_snv_GATK import ProcessGATKSNV
from garmire_SNV_calling.process_snv_GATK import PICARD_DIR
from garmire_SNV_calling.process_snv_GATK import PLATEFORM
from garmire_SNV_calling.process_snv_GATK import ORGANISM
from garmire_SNV_calling.process_snv_GATK import REF_GENOME
from garmire_SNV_calling.process_snv_GATK import DBSNP
from garmire_SNV_calling.process_snv_GATK import VCF_RESOURCES
from garmire_SNV_calling.process_snv_GATK import PROCESS_ID

from garmire_SNV_calling.config import PATH_OPOSSUM
from garmire_SNV_calling.config import PATH_FREEBAYES

from garmire_SNV_calling.config import OUTPUT_PATH_GATK


if "--do_both_callers" in argv:
    DO_BOTH_CALLERS = True
else:
    DO_BOTH_CALLERS = False

if "--path_to_data" in argv:
    PATH_TO_DATA = argv[
        argv.index("--path_to_data") + 1]
    PATH_OUTPUT =  PATH_TO_DATA + '/freebayes/'
else:
    from garmire_SNV_calling.config import PATH_OUTPUT
    from garmire_SNV_calling.config import OUTPUT_PATH_FREEBAYES


def main():
    """ """
    process_freebayes = ProcessFreebayesCaller(id=PROCESS_ID)
    if DO_BOTH_CALLERS:
        process_freebayes.process_ALL_callers()
    else:
        process_freebayes.process()


class ProcessFreebayesCaller(ProcessGATKSNV):
    """ """
    def __init__(self,
                 output_path=OUTPUT_PATH_FREEBAYES,
                 path_to_data=PATH_OUTPUT,
                 picard_dir=PICARD_DIR,
                 plateform=PLATEFORM,
                 organism=ORGANISM,
                 path_freebayes=PATH_FREEBAYES,
                 ref_genome=REF_GENOME,
                 samtools=SAMTOOLS,
                 dbsnp=DBSNP,
                 vcf_resources=VCF_RESOURCES,
                 output_path_gatk=OUTPUT_PATH_GATK,
                 respath_gatk=None,
                 **kwargs):
        """ """
        self.output_path_gatk = output_path_gatk
        self.path_freebayes =path_freebayes
        self.samtools = samtools

        ProcessGATKSNV.__init__(
            self,
            output_path=output_path,
            path_to_data=path_to_data,
            picard_dir=picard_dir,
            plateform=plateform,
            organism=organism,
            ref_genome=ref_genome,
            dbsnp=dbsnp,
            vcf_resources=vcf_resources,
            **kwargs)

        self.respath_gatk = respath_gatk

    def process(self, srr_to_process=None):
        """
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
        self._process_freebayes('recal.bam')
        self._finish_process()
        self._rm_tmp_file()

    def process_ALL_callers(self, srr_to_process=None):
        """
        """
        if srr_to_process:
            self.srr_to_process = srr_to_process

        msg = self._init_process()

        if msg:
            print(msg)
            self.stdout.write(msg)
            return

        self._init_process_gatk()
        self._launch_picard_readgroups()
        self._launch_picard_markduplicates()
        self._launch_gatk_cigar()
        self._launch_gatk_realigner_target_creator()
        self._launch_gatk_realigner_indel()
        self._launch_gatk_base_recalibrator()
        self._launch_gatk_print_reads()
        self._process_freebayes('recal.bam')
        self._launch_gatk_variant_calling()
        self._launch_gatk_variant_filtering()

        self._finish_process(ext="_GATK", out="_GATK")
        self._finish_process(ext="_freebayes", out="_freebayes")
        self._rm_tmp_file()

    def _init_process_gatk(self):
        """
        """
        if not self.respath_gatk:
            self.respath_gatk = self.output_path_gatk + \
                                "/data/" + self.srr_to_process

    def _process_samtools_calmd(self, bam_input="Aligned.sortedByCoord.out.bam"):
        """
        """
        if self.check_if_output_exists(
            "{0}/md.bam".format(self.tmppath)):
            return

        self._run_cmd(
            'echo "\n\n######## LAUNCHING SAMTOOLS CALMD ########\n"')

        cmd = "{0} calmd -b {1}/{2} {3} > {1}/md.bam".format(
            self.samtools,
            self.tmppath,
            bam_input,
            self.ref_genome)

        self._run_cmd(cmd)

        cmd = "{0} index {1}/md.bam".format(
            self.samtools,
            self.tmppath)

        self._run_cmd(cmd)

    def _process_opossum(self, bam_input="md.bam"):
        """
             " --SoftClipsExist True --KeepMismatches True " \
        """
        if self.check_if_output_exists(
            "{0}/opossum.bam".format(self.tmppath)):
            return

        self._run_cmd(
            'echo "\n\n######## LAUNCHING opossum ########\n"')

        cmd = "{0} {1}/Opossum.py --BamFile {2}/{3} " \
              " --OutFile {2}/clean.bam ".format(
            PYTHON,
            PATH_OPOSSUM,
            self.tmppath,
            bam_input
            )

        self._run_cmd(cmd)

    def _process_freebayes(self, bam_input="clean.bam"):
        """
        """
        if self.check_if_output_exists(
            "{0}/snv_filtered_freebayes.vcf".format(self.tmppath)):
            return

        self._run_cmd(
            'echo "\n\n######## LAUNCHING freebayes ########\n"')

        start_time = time()

        cmd = "{0}  -f {1} {2}/{3} > {2}/snv_filtered_freebayes.vcf".format(
            self.path_freebayes,
            self.ref_genome,
            self.tmppath,
            bam_input
            )

        self._run_cmd(cmd)

        self._run_cmd(
            'echo "\n## freebayes done in {0} s##\n"'.format(time() - start_time))

        assert(isfile("{0}/snv_filtered_freebayes.vcf".format(self.tmppath)))

if __name__ == '__main__':
    main()
