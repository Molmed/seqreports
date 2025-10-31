#!/usr/bin/env python

import pytest
import os.path
import subprocess
from bs4 import BeautifulSoup


# Run pipeline in test mode, this is done once per test session
@pytest.fixture(scope="session", autouse=True)
def result_dir(request, tmpdir_factory):
    demultiplexer = "bcl2fastq"
    if hasattr(request, "param"):
        demultiplexer = request.param

    result_dir = tmpdir_factory.mktemp("results")
    cmd = [
        "nextflow",
        "run",
        "main.nf",
        "-profile",
        "dev,test,singularity",
        "--demultiplexer",
        demultiplexer,
        "--result_dir",
        result_dir,
    ]
    if demultiplexer == "bclconvert":
        cmd = [*cmd, *["--run_folder", "test_data/230825_M04034_0043_000000000-L6NVV"]]

    subprocess.run(cmd, check=True)

    yield result_dir


def test_results_dirs_exist(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    projects_dir = os.path.join(result_dir, "projects")

    assert os.path.isdir(flowcell_dir)
    assert os.path.isdir(projects_dir)


def test_project_dirs_exist(result_dir):
    projects_dir = os.path.join(result_dir, "projects")
    projects = ["Zymo", "Qiagen", "NoProject"]

    for project in projects:
        assert os.path.isdir(os.path.join(projects_dir, project))


def test_flowcell_report_exist(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    report_path = os.path.join(
        flowcell_dir, "210510_M03910_0104_000000000-JHGJL_multiqc_report.html"
    )

    assert os.path.isfile(report_path)


def test_project_reports_exist(result_dir):
    projects_dir = os.path.join(result_dir, "projects")
    projects = ["Zymo", "Qiagen", "NoProject"]

    for project in projects:
        report_path = os.path.join(
            projects_dir,
            project,
            "210510_M03910_0104_000000000-JHGJL_" + project + "_multiqc_report.html",
        )
        assert os.path.isfile(report_path)


def check_sections_in_report(report_path, sections):
    with open(report_path, "r") as html_file:
        parser = BeautifulSoup(html_file.read(), "lxml")
        for section in sections:
            hits = parser.find_all(href="#" + section)
            assert len(hits) > 0


def test_all_sections_included_in_flowcell_report(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    report_path = os.path.join(
        flowcell_dir, "210510_M03910_0104_000000000-JHGJL_multiqc_report.html"
    )
    sections = [
        "general_stats",
        "rrna",
        "sequencing_metadata",
        "bcl2fastq",
        "interop",
        "fastq_screen",
        "fastqc",
    ]

    check_sections_in_report(report_path, sections)


@pytest.mark.parametrize("result_dir", ["bclconvert"], indirect=True)
def test_all_sections_included_in_bclcovert_flowcell_report(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    report_path = os.path.join(
        flowcell_dir, "230825_M04034_0043_000000000-L6NVV_multiqc_report.html"
    )
    sections = [
        "general_stats",
        "rrna",
        "sequencing_metadata",
        "bclconvert",
        "interop",
        "fastq_screen",
        "fastqc",
    ]

    check_sections_in_report(report_path, sections)


def test_all_sections_included_in_project_reports(result_dir):
    projects_dir = os.path.join(result_dir, "projects")
    projects = ["Zymo", "Qiagen", "NoProject"]
    sections = [
        "general_stats",
        "rrna",
        "sequencing_metadata",
        "fastq_screen",
        "fastqc",
    ]

    for project in projects:
        report_path = os.path.join(
            projects_dir,
            project,
            "210510_M03910_0104_000000000-JHGJL_" + project + "_multiqc_report.html",
        )
        check_sections_in_report(report_path, sections)
