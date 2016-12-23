# Download SRA files from NCBI (GEO)

This module provides scripts to  download and extract SRA files for High-throughput genomic data from NCBI (GEO portal) using NCBI .soft file


## Requirements
* [python 2 (>=2.7)](https://www.python.org/download/releases/2.7.2/)
* The only external software needed is [fastq-dump](http://ncbi.github.io/sra-tools/install_config.html) to extract the .sra files. Path toward the executable must be given to the config file or parsed as argument
* A folder with the name of the project must be created and the absolute path toward that folder must be given to the config file or parsed as argument
* The .soft file from [NCBI GEO](http://www.ncbi.nlm.nih.gov/geo/) website file related to the project must be downloaded (and put into the project folder (default))
  * link for dataset description GEO webpage [example](http://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85183/soft/)
  * An example soft file is also available in the ./example/ folder of the repository (default folder)

## configuration
* all global variables can be set into the file ./garmire_download_ncbi_sra/config.py or parsed as function attributes
* arguments description can be found at any time by invoking the -h (or -H) option or by consulting the config file:

```text
-PROJECT_NAME    The name of the project (defining the name of the folder)
-PATH_DATA    The absolute path where the project will be created and the SRA files downloaded and extracted
-PATH_SOFT    path toward the .soft file (with the corresponding ftp addresses for the .sra files)
-NB_THREADS    number of threads (download in parallel) to use for downloading rsa files (default 2)
-FASTQ_DUMP    path to the fastq-dump software
-FASTQ_DUMP_OPTION    options to use to extract the sra (using fastq-dump) "--split-3 -B is the default" and it is strongly recommended to keep it
-LIMIT    define the maximum number of sra files to be downloaded (default None)
```

## usage
move to folder of the git project (https://github.com/lanagarmire/SSrGE.git)

```bash
cd SSrGE
```

* Setting the global variables into the config file (download_ncbi_sra/config.py) or parsing them each time as arguments
* [optional] Running the tests:

```bash
  python ./test/test_download.py -v
  ```

* download and extract data (download by default .sra files from the example .soft file):

```bash
python garmire_download_ncbi_sra/download_data.py
```
* download and extract data (with parsing options):

```bash
python garmire_download_ncbi_sra/download_data.py -NB_THREADS 5 -PATH_SOFT tutut/...
```
* extract SRA file

```bash
python garmire_download_ncbi_sra/extract_data.py
```
* remove SRA file

```bash
python garmire_download_ncbi_sra/remove_sra.py
```

## contact and credentials
* Developer: Olivier Poirion (PhD)
* contact: opoirion@hawaii.edu