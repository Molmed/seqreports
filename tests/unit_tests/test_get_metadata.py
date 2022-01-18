#!/usr/bin/env python

import pytest
import os
import sys

# Add pipeline scripts dir to sys.path so that they can be imported
bin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin")
sys.path.append(bin_path)

from get_metadata import RunfolderInfo


@pytest.fixture
def runfolder_info():
    return RunfolderInfo("test_data/210510_M03910_0104_000000000-JHGJL", "Unaligned")    

@pytest.fixture
def nested_dictionary():
    d = {"A_key": "A_val", 
         "B_key": {"C_key": "C_val",
                   "D_key": { "E_key": "E_val"}}}
    return d

def test_find_function(runfolder_info, nested_dictionary):
    value = list(runfolder_info.find(nested_dictionary, "E_key"))
    assert value[0] == "E_val"

def test_read_run_parameters(runfolder_info):
    run_parameters = runfolder_info.read_run_parameters()
    assert len(run_parameters["RunParameters"]) == 63

def test_read_stats_json(runfolder_info):
    stats_json = runfolder_info.read_stats_json("Unaligned")
    assert len(stats_json) == 6

def test_bcl2fastq_version(runfolder_info):
    bcl2fastq_version = runfolder_info.get_bcl2fastq_version("test_data/210510_M03910_0104_000000000-JHGJL")
    assert bcl2fastq_version == "2.20.0.422"

def test_get_run_parameters(runfolder_info):
    filtered_run_parameters = runfolder_info.get_run_parameters()
    assert len(filtered_run_parameters) == 5

def test_get_read_cycles(runfolder_info):
    read_cycles = runfolder_info.get_read_cycles()
    assert len(read_cycles) == 4
    assert sum("Read" in key for key in read_cycles.keys()) == 2
    assert sum("Index" in key for key in read_cycles.keys()) == 2

def test_get_info(runfolder_info):
    results = runfolder_info.get_info()
    assert len(results) == 10
