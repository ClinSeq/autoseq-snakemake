#!/usr/bin/env python

import unittest
from unittest import mock

from pipeline.utils import utils
from pipeline.utils.clinseq_barcodes import *


class TestClinseqBarcodes(unittest.TestCase):
    def setUp(self):
        self.test_capture1 = \
            UniqueCapture("LB", "P-00000001", "CFDNA", "01234567", "TP", "CM")
        self.test_capture2 = \
            UniqueCapture("LB", "P-00000001", "CFDNA", "01234567", "TP", "WG")
        self.library = 'LB-P-NA12877-N-03098850-TD1-C31'
        self.libdir = 'tests/libraries'
        self.sampledata = { 
                            "sdid": "P-NA12877",
                            "T": "",
                            "N": ["LB-P-NA12877-N-03098850-TD1-C31"],
                            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-C31"]
                        }

    def test_extract_unique_capture_valid1(self):
        self.assertEqual(self.test_capture1,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_extract_unique_capture_valid2(self):
        self.assertEqual(self.test_capture1,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP1-CM1"))

    def test_extract_unique_capture_valid3(self):
        self.assertEqual(self.test_capture2,
                          extract_unique_capture("LB-P-00000001-CFDNA-01234567-TP1-WGS"))

    def test_extract_unique_capture_invalid(self):
        self.assertRaises(ValueError, lambda: extract_unique_capture("an_invalid_barcode"))

    def test_clinseq_barcode_is_valid_valid1(self):
        self.assertTrue(clinseq_barcode_is_valid("LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_clinseq_barcode_is_valid_invalid_too_few(self):
        self.assertFalse(clinseq_barcode_is_valid("LB-00000001-CFDNA-01234567-TP-CM2017001022000"))

    def test_clinseq_barcode_is_valid_invalid_project(self):
        self.assertFalse(clinseq_barcode_is_valid("01-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    def test_sdid_valid_valid1(self):
        self.assertTrue(sdid_valid("P-A"))

    def test_sdid_valid_invalid1(self):
        self.assertFalse(sdid_valid("P-$"))

    def test_prep_id_valid_valid1(self):
        self.assertTrue(prep_id_valid("AA1"))

    def test_prep_id_valid_invalid_too_short(self):
        self.assertFalse(prep_id_valid("AA"))

    def test_prep_id_valid_invalid_no_number(self):
        self.assertFalse(prep_id_valid("AAZ"))

    def test_capture_id_valid_valid_wgs(self):
        self.assertTrue(capture_id_valid("WGS"))

    def test_capture_id_valid_invalid_no_number(self):
        self.assertFalse(capture_id_valid("AAZ"))

    def test_capture_id_valid_invalid_too_short(self):
        self.assertFalse(capture_id_valid("AA"))

    def test_extract_kit_id_valid(self):
        self.assertEqual(extract_kit_id("SomeString"), "So")

    def test_extract_kit_id_invalid(self):
        self.assertRaises(ValueError, lambda: extract_kit_id("C"))

    def test_extract_clinseq_barcodes_invalid_file_extension(self):
        self.assertRaises(ValueError, lambda: extract_clinseq_barcodes("some_file.invalid_extension"))

    def test_validate_clinseq_barcodes_two_valid(self):
        validate_clinseq_barcodes(["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
                                   "LB-P-00000002-T-01234567-TP201701011541-CM2017001022001"])

    def test_validate_clinseq_barcodes_none(self):
        validate_clinseq_barcodes([])

    def test_validate_clinseq_barcodes_one_invalid(self):
        self.assertRaises(ValueError, lambda: validate_clinseq_barcodes(["an_invalid_barcode"]))


    def test_validate_clinseq_barcodes_one_of_two_invalid(self):
        self.assertRaises(ValueError, lambda: validate_clinseq_barcodes([
            "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
            "an_invalid_barcode"]))

    def test_convert_barcodes_to_sampledict_invalid_barcode(self):
        self.assertRaises(ValueError, lambda: convert_barcodes_to_sampledict([
            "an_invalid_barcode", "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"]))

    def test_convert_barcodes_to_sampledict_two_barcodes_same_sdid(self):
        sample_dict = convert_barcodes_to_sampledict([
            "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
            "LB-P-00000001-N-01234568-TP201701011540-CM2017001022000"])
        self.assertEqual(list(sample_dict.keys()), ["P-00000001"])

    def test_convert_barcodes_to_sampledict_two_barcodes_different_sdid(self):
        sample_dict = convert_barcodes_to_sampledict([
            "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
            "LB-P-00000002-N-01234568-TP201701011540-CM2017001022000"])
        self.assertEqual(set(sample_dict.keys()), set(["P-00000001", "P-00000002"]))
    
    def test_extract_clinseq_barcodes_txt(self):
        mocked_open = mock.mock_open(read_data='a_mock_barcode\nanother_mock_barcode\n')
        with mock.patch('pipeline.utils.clinseq_barcodes.open', mocked_open, create=True):
            observed_len = len(extract_clinseq_barcodes("test.txt"))
            self.assertEqual(observed_len, 2)

    def test_find_fastqs_for_no_library(self):
        """
        test that find_fastqs return (None,None) if called with library=None
        """
        files = find_fastqs(library=None, libdir=self.libdir)
        self.assertEqual(files, (None, None))

    def test_find_fastq_gz(self):
        """
        test that files on the format *_1.fastq.gz / *_2.fastq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('2_LB-P-NA12877-N-03098850-TD1-TT1_1.fastq.gz', files_basenames)
        self.assertIn('2_LB-P-NA12877-N-03098850-TD1-TT1_2.fastq.gz', files_basenames)
    
    @mock.patch('pipeline.utils.clinseq_barcodes.os.path.exists')
    @mock.patch('pipeline.utils.clinseq_barcodes.find_fastqs')
    def test_data_available_for_clinseq_barcode_file_exists(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = True
        mock_os_path_exists.return_value = True
        self.assertTrue(
            data_available_for_clinseq_barcode("test_libdir",
                                               "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))
    
    @mock.patch('pipeline.utils.clinseq_barcodes.os.path.exists')
    @mock.patch('pipeline.utils.clinseq_barcodes.find_fastqs')
    def test_data_available_for_clinseq_barcode_no_dir(self, mock_find_fastqs, mock_os_path_exists):
        mock_find_fastqs.return_value = True
        mock_os_path_exists.return_value = False
        self.assertFalse(
            data_available_for_clinseq_barcode("test_libdir",
                                               "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"))

    @mock.patch('pipeline.utils.clinseq_barcodes.data_available_for_clinseq_barcode')
    def test_check_sampledata(self, mock_data_available_for_clinseq_barcode):
        """
        test that all clinseq barcodes are valid and fq exists
        """
        mock_data_available_for_clinseq_barcode.return_value = True
        barcodes = ["LB-P-NA12877-N-03098850-TD1-C31", "LB-P-NA12877-CFDNA-03098850-TD1-C31"]
        samples, clinseq_barcodes = check_sampledata(self.libdir, self.sampledata)
        self.assertEqual(self.sampledata, samples)
        self.assertEqual(barcodes, clinseq_barcodes)

    @mock.patch('pipeline.utils.clinseq_barcodes.data_available_for_clinseq_barcode')
    def test_check_sampledata_not_valid(self, mock_data_available_for_clinseq_barcode):
        """
        test that all clinseq barcodes are valid and fq exists
        """
        mock_data_available_for_clinseq_barcode.return_value = False
        barcodes = ["LB-P-NA12877-N-03098850-TD1-C31", "LB-P-NA12877-CFDNA-03098850-TD1-C31"]
        self.assertRaises(FileNotFoundError, lambda: check_sampledata(self.libdir, self.sampledata))
    
    