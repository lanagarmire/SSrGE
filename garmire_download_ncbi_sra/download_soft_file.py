from garmire_SNV_calling.bash_utils import exec_cmd

from garmire_download_ncbi_sra.config import PATH_DATA
from garmire_download_ncbi_sra.config import PROJECT_NAME
from garmire_download_ncbi_sra.config import SOFT_ID

from glob import glob

from collections import Counter

from os import mkdir

from os.path import isfile

import re

from collections import defaultdict

from os.path import isdir

from datetime import datetime
import json


def main():
    """ """
    if not isdir(PATH_DATA):
        mkdir(PATH_DATA)

    download_and_process_soft(SOFT_ID)


def download_and_process_soft(gse, erase=False):
    """
    """
    if not erase:
        if glob('{0}/{1}*'.format(PATH_DATA, gse)):
            print('soft file seems existing for: {0}'.format(gse))
            return

    print('downloadin: {0}...'.format(gse))

    address = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{0}nnn/{1}/soft/{1}_family.soft.gz'.format(
        gse[:-3],
        gse
    )

    exec_cmd('wget {0} -O {1}/{2}.soft.gz'.format(address, PATH_DATA, PROJECT_NAME))
    exec_cmd('gzip -d {0}/{1}.soft.gz'.format(PATH_DATA, PROJECT_NAME))

    read_soft('{0}/{1}.soft'.format(PATH_DATA, PROJECT_NAME))

def read_soft(soft_file):
    """
    """
    gse_dict = extract_gsm_from_soft(soft_file)

    if not gse_dict:
        print('soft_file:{0} empty!'.format(soft_file))
        return

    n_samples = len(gse_dict)

    organism = Counter([gse_dict[gse]['organism_code'] for gse in gse_dict])

    organism = sorted(organism.items(), key=lambda x:x[1], reverse=True)[0][0]
    organism = organism.split()[0]

    f_stat = open('{0}/statistics.json'.format(PATH_DATA), 'w')
    f_meta = open('{0}/metadata.json'.format(PATH_DATA), 'w')

    f_meta.write(json.dumps(gse_dict, indent=2))

    f_stat.write(json.dumps({
        'organism':organism,
        "nb_samples": n_samples
    }, indent=2))

    print("organism: {0}".format(organism))
    print("number of samples: {0}".format(n_samples))


def extract_gsm_from_soft(
        soft_file,
        flatten_gsm=False,
        remove_not_sra=True):
    """
    """
    gse_dict = {}

    assert(isfile(soft_file))

    f_soft = open(soft_file)
    line = f_soft.readline()

    while line:
        if line.count('^SAMPLE'):
            data = defaultdict(list)
            gse = line.strip('\n').split(' = ', 1)[1].strip()
            data['GSE'] = gse

            line = f_soft.readline()

            while line and line[0] == '!':
                key, value = line.split(' = ', 1)

                key = key.strip('! ')
                key = key.replace('/', '_')
                value = value.strip('\n ')

                if key[:7] == 'Sample_':
                    key = key[7:]

                if value[:6] == 'ftp://':
                    data['ftp'].append(value)

                sra = re.findall('https://www.ncbi.nlm.nih.gov/sra?term=SRX[0-9]+', value)

                geo_organism = data['organism_ch1']

                if geo_organism:
                    geo_organism = geo_organism[0]

                    if geo_organism == 'Homo sapiens' or geo_organism == 'Homo':
                        data['organism_code'] = 'HUMAN'
                    elif geo_organism == 'Mus musculus' or geo_organism == 'Mus':
                        data['organism_code'] = 'MOUSE'

                if sra:
                    data['SRA'].append(sra[0])

                data[key].append(value)

                line = f_soft.readline()

            if 'relation' in data:
                for relation in data['relation']:
                    key, value = relation.split(': ')
                    key.strip(), value.strip()
                    data[key] = value

            if flatten_gsm:
                for key in data:
                    data[key] = check_value(key, data[key])

                    if len(data[key]) == 1:
                        data[key] = data[key][0]

            if 'SRA' in data or not remove_not_sra:
                gse_dict[gse] = data

        else:
            line = f_soft.readline()

    return gse_dict


def check_value(key, values):
    """
    """
    is_list= True

    if not isinstance(values, list):
        is_list = False
        values = [values]

    values = map(format_value, values)

    if key.count('zip') and not isinstance(values[0], int):
        values = map(lambda x:0, values)

    if key.count('phone') and not isinstance(values[0], int):
        values = map(lambda x:0, values)

    if not is_list:
        values = values[0]

    return values

def format_value(value):
    """
    """
    if value.isdigit():
        value = int(value)
    elif re.findall('[A-Z][a-z][a-z] [0-9]{2} [0-9]{4}', value):
        value = re.findall('[A-Z][a-z][a-z] [0-9]{2} [0-9]{4}', value)[0]
        value = datetime.strptime(value, '%b %d %Y')

    return value


if __name__ == '__main__':
    main()
