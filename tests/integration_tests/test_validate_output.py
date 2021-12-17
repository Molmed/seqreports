#!/usr/bin/env python

import pytest
import os.path
import subprocess
from bs4 import BeautifulSoup
import itertools


@pytest.fixture(scope="session", autouse=True)
def result_dir(tmpdir_factory):
    result_dir = tmpdir_factory.mktemp("results")
    subprocess.run(
        [
            "nextflow",
            "run",
            "main.nf",
            "-profile",
            "dev,test,singularity",
            "--result_dir",
            result_dir,
        ],
        check=True,
    )
    yield result_dir


@pytest.fixture
def flowcell_report_dir(result_dir):
    return os.path.join(result_dir, "flowcell_report")


@pytest.fixture
def project_reports_dir(result_dir):
    return os.path.join(result_dir, "projects")


@pytest.fixture
def projects():
    return ["Zymo", "Qiagen", "NoProject"]


@pytest.fixture
def flowcell_report(flowcell_report_dir):
    return os.path.join(
        flowcell_report_dir, "210510_M03910_0104_000000000-JHGJL_multiqc_report.html"
    )


@pytest.fixture
def project_reports(project_reports_dir, projects):
    report_list = []
    for project in projects:
        report_list.append(
            os.path.join(
                project_reports_dir,
                project,
                "210510_M03910_0104_000000000-JHGJL_"
                + project
                + "_multiqc_report.html",
            )
        )
    return report_list


@pytest.fixture
def flowcell_report_sections():
    sections = [
        "general_stats",
        "rrna",
        "sequencing_metadata",
        "bcl2fastq",
        "interop",
        "fastq_screen",
        "fastqc",
    ]
    return sections


@pytest.fixture
def project_report_sections():
    sections = [
        "general_stats",
        "rrna",
        "sequencing_metadata",
        "fastq_screen",
        "fastqc",
    ]
    return sections


def test_results_dirs_exist(flowcell_report_dir, project_reports_dir):
    assert os.path.isdir(flowcell_report_dir)
    assert os.path.isdir(project_reports_dir)


def test_project_dirs_exist(project_reports_dir, projects):
    for project in projects:
        assert os.path.isdir(os.path.join(project_reports_dir, project))


def test_reports_exist(flowcell_report, project_reports):
    reports = project_reports + [flowcell_report]
    for report in reports:
        assert os.path.isfile(report)


def test_all_sections_included(
    flowcell_report,flowcell_report_sections, project_reports, project_report_sections
):
    def check_sections_in_reports(reports, sections):
        for report_path in reports:
            with open(report_path, "r") as html_file:
                parser = BeautifulSoup(html_file.read(), "lxml")
                for section in sections:
                    hits = parser.find_all(href="#" + section)
                    assert len(hits) > 0

    check_sections_in_reports([flowcell_report], flowcell_report_sections)
    check_sections_in_reports(project_reports, project_report_sections)
