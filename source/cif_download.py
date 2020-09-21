import urllib.request
import os
from urllib.error import HTTPError
from pathlib import Path
from .zeotype import Zeotype

def download_cif(code: str, data_dir='data'):
    Path(data_dir).mkdir(parents=True, exist_ok=True)  # create output directory if it doesn't exist
    output_path = os.path.join(data_dir, code + '.cif')  # sets ouput path to data_dir/code.cif
    root_url = "https://europe.iza-structure.org/IZA-SC/cif/"  # sets root URL for the
    url = root_url + code + '.cif'
    try:
        urllib.request.urlretrieve(url, output_path)
    except HTTPError as err:
        if err.code == 404:
            print("error code 404: Specified Zeolite Framework not found")
            raise
        else:
            raise


def build_zeolite_from_code(code, data_dir='data'):
    output_path = os.path.join(data_dir, code + '.cif')
    if not os.path.exists(output_path):
        download_cif(code, data_dir)  # will throw error if not successful
    return Zeotype.build_zeolite_from_cif(output_path)