#!/usr/bin/env python

import unittest

from unittest.mock import patch, mock_open
from pipeline.utils import utils
import snakemake

class TestWorkflow(unittest.TestCase):

    def setUp(self):
        self.snakefile = "pipeline/autoseq/Snakefile"
        self.config = "config.yml"
        self.reference = "tests/dummy_genome/dummy_genome.json"
    
    @patch("pipeline.utils.utils.os.path.isfile")
    @patch("pipeline.utils.utils.get_chromosomes")
    def test_autoseq_valid(self, mock_isfile, mock_get_chromosomes):
        mock_isfile.return_value = True
        mocked_open = mock_open(read_data='{"1", "2", "3", "4", "5", "6", "7", "8", "X", "Y"}')
        mock_get_chromosomes.return_value = {'1', '2', '3', '4', '5', '6', '7', '8', 'X', 'Y'}

        with patch("pipeline.utils.utils.open", mocked_open, create=True):
            self.assertTrue(snakemake.snakemake(self.snakefile,
                                                configfiles=[self.config],
                                                dryrun=True))

    @patch("pipeline.utils.utils.os.path.isfile")
    @patch("pipeline.utils.utils.get_chromosomes")
    def test_autoseq_invalid(self, mock_isfile, mock_get_chromosomes):
        mock_isfile.return_value = True
        mocked_open = mock_open(read_data=None)
        mock_get_chromosomes.return_value = None

        with patch("pipeline.utils.utils.open", mocked_open, create=True):
            self.assertFalse(snakemake.snakemake(self.snakefile,
                                                configfiles=[self.config],
                                                dryrun=True))
