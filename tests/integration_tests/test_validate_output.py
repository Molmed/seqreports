#!/usr/bin/env python

import pytest
import os.path
import subprocess
import tempfile


@pytest.fixture(scope="session", autouse=True)
def result_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        result_dir = os.path.join(tmpdir, "results")
        test_run = subprocess.run(
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


def test_results_dirs_exist(flowcell_report_dir, project_reports_dir):
    assert os.path.isdir(flowcell_report_dir)
    assert os.path.isdir(project_reports_dir)


def test_project_dirs_exist(project_reports_dir, projects):
    for project in projects:
        assert os.path.isdir(os.path.join(project_reports_dir, project))


def test_reports_exist(flowcell_report_dir, project_reports_dir):
    assert os.path.isfile(
        os.path.join(
            flowcell_report_dir,
            "210510_M03910_0104_000000000-JHGJL_multiqc_report.html",
        )
    )
