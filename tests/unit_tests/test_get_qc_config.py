#!/usr/bin/env python

import pytest
import os
import sys

# Add pipeline scripts dir to sys.path so that they can be imported
bin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin")
sys.path.append(bin_path)

from get_qc_config import convert_to_dict, convert_to_multiqc_config


@pytest.fixture
def checkqc_config():
    return [
        {"name": "ClusterPFHandler", "warning": 18, "error": "unknown"},
        {"name": "Q30Handler", "warning": 80, "error": "unknown"},
        {
            "name": "ErrorRateHandler",
            "allow_missing_error_rate": False,
            "warning": 2,
            "error": "unknown",
        },
        {"name": "ReadsPerSampleHandler", "warning": "unknown", "error": 13.5},
        {"name": "UndeterminedPercentageHandler", "warning": "unknown", "error": 9},
        {
            "name": "UnidentifiedIndexHandler",
            "significance_threshold": 1,
            "white_listed_indexes": [".*N.*", "G{6,}"],
        },
    ]


@pytest.fixture
def checkqc_config_dict():
    return {
        "ClusterPFHandler": {
            "name": "ClusterPFHandler",
            "warning": 18,
            "error": "unknown",
        },
        "Q30Handler": {"name": "Q30Handler", "warning": 80, "error": "unknown"},
        "ErrorRateHandler": {
            "name": "ErrorRateHandler",
            "allow_missing_error_rate": False,
            "warning": 2,
            "error": "unknown",
        },
        "ReadsPerSampleHandler": {
            "name": "ReadsPerSampleHandler",
            "warning": "unknown",
            "error": 13.5,
        },
        "UndeterminedPercentageHandler": {
            "name": "UndeterminedPercentageHandler",
            "warning": "unknown",
            "error": 9,
        },
        "UnidentifiedIndexHandler": {
            "name": "UnidentifiedIndexHandler",
            "significance_threshold": 1,
            "white_listed_indexes": [".*N.*", "G{6,}"],
        },
    }


def test_convert_to_dict(checkqc_config):
    checkqc_config_dict = convert_to_dict(checkqc_config)
    assert isinstance(checkqc_config_dict, dict)
    assert len(checkqc_config_dict.keys()) == 6


def test_convert_to_multiqc_config(checkqc_config_dict):
    multiqc_config = convert_to_multiqc_config(checkqc_config_dict)
    multiqc_section = multiqc_config["table_cond_formatting_rules"]
    assert multiqc_section["Error"]["warn"] == [{"gt": 2}]
    assert multiqc_section["mqc-generalstats-bcl2fastq-total"]["fail"] == [{"lt": 13.5}]
    assert multiqc_section["percent_Q30"]["warn"] == [{"lt": 80}]
    assert multiqc_section["total"]["warn"] == [{"lt": 18}]
