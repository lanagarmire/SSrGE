"""generate STAR GENOME INDEX"""

from sys import stdout as sys_stdout
from os import popen

from os.path import isdir
from distutils.dir_util import mkpath

from garmire_SNV_calling.config import PATH_STAR_SOFTWARE
from garmire_SNV_calling.config import STAR_INDEX_PATH
from garmire_SNV_calling.config import ANNOTATION_PATH
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import STAR_THREADS
from garmire_SNV_calling.config import STAR_INDEX_READ_LENGTH


def main():
    """ """
    star_index_path = "{0}READ{1}/".format(STAR_INDEX_PATH.rstrip('/'),
                                          STAR_INDEX_READ_LENGTH)
    print "######## computing STAR index ########\npath:{0}\n"\
        .format(star_index_path)

    if not isdir(star_index_path):
        mkpath(star_index_path)

    cmd = "{0} --runMode genomeGenerate --runThreadN {1}"\
          " --genomeDir {2} --genomeFastaFiles {3} --sjdbGTFfile {4}"\
          " --sjdbOverhang {5}"\
          .format(
              PATH_STAR_SOFTWARE,
              STAR_THREADS,
              star_index_path,
              REF_GENOME,
              ANNOTATION_PATH,
              STAR_INDEX_READ_LENGTH
    )
    stdout = popen(cmd)
    c = stdout.read(1)

    while c:
        sys_stdout.write(c)
        sys_stdout.flush()
        c = stdout.read(1)

if __name__ == "__main__":
    main()
