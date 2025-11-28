#!/usr/bin/env python

import pytest
import os.path
import subprocess
from bs4 import BeautifulSoup


# Run pipeline in test mode, this is done once per test session
@pytest.fixture(scope="session", autouse=True)
def result_dir(request, tmpdir_factory):
    demultiplexer = request.param

    result_dir = tmpdir_factory.mktemp("results")
    extra_profile = "test_bclconvert" if demultiplexer == "bclconvert" else "test"

    subprocess.run(
        [
            "nextflow",
            "run",
            "main.nf",
            "-profile",
            f"dev,{extra_profile},singularity",
            "--demultiplexer",
            demultiplexer,
            "--result_dir",
            result_dir,
        ],
        check=True,
    )

    yield result_dir


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
def test_results_dirs_exist(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    projects_dir = os.path.join(result_dir, "projects")

    assert os.path.isdir(flowcell_dir)
    assert os.path.isdir(projects_dir)


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
def test_project_dirs_exist(result_dir):
    projects_dir = os.path.join(result_dir, "projects")
    projects = ["Zymo", "Qiagen", "NoProject"]

    for project in projects:
        assert os.path.isdir(os.path.join(projects_dir, project))


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
def test_flowcell_report_exist(result_dir):
    flowcell_dir = os.path.join(result_dir, "flowcell_report")
    report_path = os.path.join(
        flowcell_dir, "210510_M03910_0104_000000000-JHGJL_multiqc_report.html"
    )

    assert os.path.isfile(report_path)


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
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


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
def check_sections_in_report(report_path, sections):
    with open(report_path, "r") as html_file:
        parser = BeautifulSoup(html_file.read(), "lxml")
        for section in sections:
            hits = parser.find_all(href="#" + section)
            assert len(hits) > 0


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
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


@pytest.mark.parametrize("result_dir", ["bcl2fastq"], indirect=True)
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
