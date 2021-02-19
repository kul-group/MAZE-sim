from unittest import TestCase
from maze.cif_download import download_cif
import os
from pathlib import Path

class test_cif_download(TestCase):
    def test_download_cif(self):
        with self.subTest(msg='CHA cif download default arguments'):
            expected_output_file = os.path.join('data', 'CHA.cif')
            if os.path.exists(expected_output_file):  # if file exists, delete it
                os.remove(expected_output_file)
            download_cif('CHA')
            self.assertTrue(os.path.exists(expected_output_file))
        with self.subTest(msg='CHA cif download specified folder location'):
            expected_output_dir = os.path.join('data','CHA_test')
            expected_output_file = os.path.join(expected_output_dir,'CHA.cif')
            Path(expected_output_file).parents[0].unlink(missing_ok=True)  # delete folder if it exists
            download_cif('CHA', data_dir=expected_output_dir)
            self.assertTrue(os.path.exists(expected_output_file))

