
import unittest
from unittest import mock
from pipeline.utils.utils import *
from pipeline.utils.clinseq_barcodes import UniqueCapture


class TestUtils(unittest.TestCase):
    def setUp(self):
        self.test_capture = \
            UniqueCapture("LB", "P-NA12877", "N", "03098850", "TD1", "C31")
        self.bamfiles = ["LB-P-NA12877-N-03098850-TD1-C31_nodups.bam"]

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
    
    def test_get_scheduler(self):
        script_path = get_scheduler('slurm', 'pyscript')
        config_path = get_scheduler('slurm', 'config')
        self.assertTrue(script_path.endswith("slurm_submit.py"))
        self.assertTrue(config_path.endswith("cluster_config.json"))
    
    def test_get_scheduler_invalid(self):
        self.assertRaises(FileNotFoundError, lambda: get_scheduler('not_valid', 'pyscript'))
    
    def test_get_capture_bam(self):
        self.assertEqual(get_capture_bam(self.test_capture, self.bamfiles), 
                        "LB-P-NA12877-N-03098850-TD1-C31_nodups.bam")
    
    def test_get_capture_bam_invalid(self):
        self.assertFalse(get_capture_bam(self.test_capture, []))


