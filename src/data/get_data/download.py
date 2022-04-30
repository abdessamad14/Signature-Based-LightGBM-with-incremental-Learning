from definitions import *
import os
from io import BytesIO
import zipfile
from urllib.request import urlopen


def main():
    """Downloads the raw sepsis data.

    There are two training sets, A and B. This function will download and unzip both sets into `data/raw/training_setA`
    and `data/raw/training_setB`.
    """

    base_dir = DATA_DIR + '/raw'
    assert os.path.isdir(DATA_DIR), 'Please make a directory at {ROOT}/data. Note: we leave this down to the user' \
                                    'to give you the option of making the directory a symlink. The estimated space ' \
                                    'required is around 3GB.'

    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)

    # Two save dirs
    save_loc_a = base_dir + '/training_setA'
    save_loc_b = base_dir + '/training_setB'
    url_a = 'https://archive.physionet.org/users/shared/challenge-2019/training_setA.zip'
    url_b = 'https://archive.physionet.org/users/shared/challenge-2019/training_setB.zip'

    locs = [(save_loc_a, url_a), (save_loc_b, url_b)]
    print('Begin file download...')
    for save_loc, url in locs:
        if os.path.exists(save_loc):
            continue

        if not os.path.exists(save_loc):
            r = urlopen(url)
            z = zipfile.ZipFile(BytesIO(r.read()))
            z.extractall(base_dir)

    if os.path.isdir(base_dir + '/training'):
        os.rename(base_dir + '/training', base_dir + '/training_setA')
    print('End file download.')

if __name__ == '__main__':
    main()
