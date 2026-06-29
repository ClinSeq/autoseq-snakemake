
import unittest
from unittest import mock
from pipeline.utils.utils import *
from pipeline.utils.clinseq_barcodes import UniqueCapture


class Wildcards:
    def __init__(self, sample):
        self.sample = sample


class TestUtils(unittest.TestCase):
    def setUp(self):
        self.test_capture = \
            UniqueCapture("LB", "P-NA12877", "N", "03098850", "TD1", "C31")
        self.bamfiles = ["LB-P-NA12877-N-03098850-TD1-C31_nodups.bam"]
        self.reference = {
            "targets": {  
                "probio_comprehensive3": {
                    "blacklist-bed": None, 
                    "cnvkit-ref": {
                        "KAPA_HYPERPREP": {
                            "CFDNA": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.CFDNA.cnn", 
                            "N": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.N.cnn", 
                            "T": "intervals/targets/probio_comprehensive3.KAPA_HYPERPREP.T.cnn"
                        }
                    }, 
                    "msings-baseline": "intervals/targets/probio_comprehensive3.msings.baseline", 
                    "msings-bed": "intervals/targets/probio_comprehensive3.msings.bed", 
                    "msings-msi_intervals": "intervals/targets/probio_comprehensive3.msings.msi_intervals", 
                    "msisites": "intervals/targets/probio_comprehensive3.slopped20.msisites.tsv", 
                    "purecn_targets": None, 
                    "targets-bed-slopped20": "intervals/targets/probio_comprehensive3.slopped20.bed", 
                    "targets-bed-slopped20-gz": "intervals/targets/probio_comprehensive3.slopped20.bed.gz", 
                    "targets-interval_list": "intervals/targets/probio_comprehensive3.interval_list", 
                    "targets-interval_list-slopped20": "intervals/targets/probio_comprehensive3.slopped20.interval_list"
                }
            }
        }

    def test_convert_to_absolute_path_non_string1(self):
        self.assertEqual(
            convert_to_absolute_path(None, "/dummy/base/dir"),
            None
        )

    def test_convert_to_absolute_path_non_string2(self):
        self.assertEqual(
            convert_to_absolute_path(1, "/dummy/base/dir"),
            1
        )

    def test_convert_to_absolute_path_absolute_path(self):
        self.assertEqual(
            convert_to_absolute_path("/an/absolute/path", "/dummy/base/dir"),
            "/an/absolute/path"
        )

    @mock.patch('pipeline.utils.utils.os.path.isfile')
    def test_convert_to_absolute_path_relative_path_file(self, mock_isfile):
        mock_isfile.return_value = True
        self.assertEqual(
            convert_to_absolute_path("a_terminal_filename", "/dummy/base/dir"),
            "/dummy/base/dir/a_terminal_filename"
        )

    @mock.patch('pipeline.utils.utils.os.path.isdir')
    def test_convert_to_absolute_path_relative_path_dir(self, mock_isdir):
        mock_isdir.return_value = True
        self.assertEqual(
            convert_to_absolute_path("a_terminal_dirname", "/dummy/base/dir"),
            "/dummy/base/dir/a_terminal_dirname"
        )

    @mock.patch('pipeline.utils.utils.os.path.isfile')
    @mock.patch('pipeline.utils.utils.os.path.isdir')
    def test_convert_to_absolute_path_relative_not_there(self, mock_isdir, mock_isfile):
        mock_isdir.return_value = False
        mock_isfile.return_value = False
        self.assertEqual(
            convert_to_absolute_path("a_relative_filename_not_existing", "/dummy/base/dir"),
            "a_relative_filename_not_existing"
        )

    @mock.patch('pipeline.utils.utils.os.path.isfile')
    def test_make_paths_absolute(self, mock_isfile):
        mock_isfile.return_value = True
        input_dict = {
            "some_key1": None,
            "some_key2": 1,
            "some_key3": "a_terminal_filename",
            "nested_dict": {"some_key4": "filename2"}
        }
        output_dict = make_paths_absolute(input_dict, "/dummy/base/dir")
        self.assertEqual(output_dict["some_key1"], None)
        self.assertEqual(output_dict["some_key2"], 1)
        self.assertEqual(output_dict["some_key3"], "/dummy/base/dir/a_terminal_filename")
        self.assertEqual(output_dict["nested_dict"]["some_key4"], "/dummy/base/dir/filename2")

    @mock.patch('pipeline.utils.utils.os.path.isfile')
    def test_make_paths_absolute(self, mock_isfile):
        mock_isfile.return_value = False
        input_dict = {
            "some_key1": None,
            "some_key2": 1,
            "some_key3": "a_terminal_filename",
            "nested_dict": {"some_key4": "filename2"}
        }
        output_dict = make_paths_absolute(input_dict, "/dummy/base/dir")
        self.assertEqual(output_dict["some_key1"], None)
        self.assertEqual(output_dict["some_key2"], 1)
        self.assertEqual(output_dict["some_key3"], "a_terminal_filename")
        self.assertEqual(output_dict["nested_dict"]["some_key4"], "filename2")
    
    def test_get_capture_bam(self):
        self.assertEqual(get_capture_bam(self.test_capture, self.bamfiles), 
                        "LB-P-NA12877-N-03098850-TD1-C31_nodups.bam")
    
    def test_get_capture_bam_invalid(self):
        self.assertFalse(get_capture_bam(self.test_capture, []))

    def test_get_capture_svs(self):
        wildcards = Wildcards("LB-P-NA12877-CFDNA-03098850-TD1-C31")
        outdir = "/dummy"
        self.assertEqual(4, len(get_capture_svs(wildcards, outdir)))
        self.assertTrue(isinstance(get_capture_svs(wildcards, outdir), dict))
        self.assertTrue(get_capture_svs(wildcards, outdir)['DEL'].endswith(".gtf"))
    
    def test_get_readgroup(self):
        wildcards = Wildcards("LB-P-NA12877-CFDNA-03098850-TD1-C31")
        readgp = '"@RG\\tID:LB-P-NA12877-CFDNA-03098850-TD1-C31\\tSM:LB-P-NA12877-CFDNA-03098850\\tLB:TD1\\tPL:ILLUMINA"'
        self.assertTrue(isinstance(get_readgroup(wildcards), str))
        self.assertEqual(get_readgroup(wildcards), readgp)
    
    def test_get_targets(self):
        wildcards = Wildcards("LB-P-NA12877-CFDNA-03098850-TD1-C31")
        self.assertEqual(get_targets(wildcards, self.reference, "targets-bed-slopped20"), 
                        "intervals/targets/probio_comprehensive3.slopped20.bed")

    def test_get_target_name(self):
        wildcards = Wildcards("LB-P-NA12877-CFDNA-03098850-TD1-C31")
        self.assertEqual(get_target_name(wildcards), 
                        "probio_comprehensive3")