#!/usr/bin/env python

import unittest

from unittest.mock import patch, mock_open
from pipeline.utils import utils
import snakemake

class TestWorkflow(unittest.TestCase):

    def setUp(self):
        self.snakefile = "pipeline/autoseq/Snakefile"
        self.to_snakefile = "pipeline/tumor_only/Snakefile"
        self.snakefile_sd = "pipeline/autoseq-sd/Snakefile"
        self.wgs_snakefile = "pipeline/autoseq-wgs/Snakefile"
        self.config = "tests/test_config.yml"
        self.config_sd = "tests/test_config_sd.yml"
        self.config_noumi = "tests/test_config_non_umi.yml"
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
    def test_autoseq_wo_umi(self, mock_isfile, mock_get_chromosomes):
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

    @patch("pipeline.utils.utils.os.path.isfile")
    @patch("pipeline.utils.utils.get_chromosomes")
    @patch("os.symlink")
    @patch("os.makedirs")
    def test_tumor_only_valid(self, mock_isfile, mock_get_chromosomes, mock_os_symlink, mock_makedirs):
        mock_isfile.return_value = True
        mocked_open = mock_open(read_data='{"1", "2", "3", "4", "5", "6", "7", "8", "X", "Y"}')
        mock_get_chromosomes.return_value = {'1', '2', '3', '4', '5', '6', '7', '8', 'X', 'Y'}
        mock_os_symlink.return_value = True
        mock_makedirs.return_value = True

        with patch("pipeline.utils.utils.open", mocked_open, create=True):
            self.assertTrue(snakemake.snakemake(self.to_snakefile,
                                                configfiles=[self.config],
                                                dryrun=True))
    
    @patch("pipeline.utils.utils.os.path.isfile")
    @patch("pipeline.utils.utils.get_chromosomes")
    def test_autoseq_sd_valid(self, mock_isfile, mock_get_chromosomes):
        mock_isfile.return_value = True
        mocked_open = mock_open(read_data='{"1", "2", "3", "4", "5", "6", "7", "8", "X", "Y"}')
        mock_get_chromosomes.return_value = {'1', '2', '3', '4', '5', '6', '7', '8', 'X', 'Y'}

        with patch("pipeline.utils.utils.open", mocked_open, create=True):
            self.assertTrue(snakemake.snakemake(self.snakefile_sd,
                                                configfiles=[self.config_sd],
                                                dryrun=True))

    @patch("pipeline.utils.utils.os.path.isfile")
    def test_autoseq_wgs_valid(self, mock_isfile):
        mock_isfile.return_value = True
        with patch("pipeline.utils.utils.open", create=True):
            self.assertTrue(snakemake.snakemake(self.wgs_snakefile,
                                                configfiles=[self.config],
                                                dryrun=True))