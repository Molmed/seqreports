#!/usr/bin/env python
from __future__ import print_function
import xmltodict
from collections import OrderedDict
import re
import glob
import csv
import argparse
import os
import json
from operator import itemgetter
from pathlib import Path
import yaml


class RunfolderInfo:
    def __init__(self, runfolder):
        self.runfolder = runfolder
        self.run_info = self.read_run_info()
        self.run_parameters = self.read_run_parameters()
        self.description_and_identifier = OrderedDict()
        self.run_parameters_tags = {
            "RunId": "Run ID",
            "RunID": "Run ID",
            "InstrumentType": "Instrument type",
            "ApplicationName": "Control software",
            "Application": "Control software",
            "ApplicationVersion": "Control software version",
            "SystemSuiteVersion": "Control software version",
            "Flowcell": "Flowcell type",
            "FlowCellMode": "Flowcell type",
            "ReagentKitVersion": "Reagent kit version",
            "RTAVersion": "RTA Version",
            "RtaVersion": "RTA Version",
        }

    def find(self, d, tag):
        if isinstance(d, dict):
            if tag in d:
                yield d[tag]
            for k, v in d.items():
                if isinstance(v, dict):
                    yield from self.find(v, tag)
                if isinstance(v, list):
                    for i in v:
                        yield from self.find(i, tag)

    def read_run_info(self):
        run_info = os.path.join(self.runfolder, "RunInfo.xml")
        if os.path.exists(run_info):
            with open(run_info) as f:
                return xmltodict.parse(f.read())
        else:
            return None

    def read_run_parameters(self):
        alt_1 = os.path.join(self.runfolder, "runParameters.xml")
        alt_2 = os.path.join(self.runfolder, "RunParameters.xml")
        if os.path.exists(alt_1):
            with open(alt_1) as f:
                return xmltodict.parse(f.read())
        elif os.path.exists(alt_2):
            with open(alt_2) as f:
                return xmltodict.parse(f.read())
        else:
            return None

    def find_flowcell_type_novaseqx(self):
        try:
            consumables = self.run_parameters["RunParameters"]["ConsumableInfo"][
                "ConsumableInfo"
            ]
            flowcell = next(
                consumable
                for consumable in consumables
                if consumable["Type"] == "FlowCell"
            )
            flowcell_type = flowcell.get("Name", flowcell.get("Mode"))
        except (KeyError, StopIteration):
            return None
        return {"Flowcell type": flowcell_type}

    def get_software_version(self, runfolder):
        pipeline_dir = Path(runfolder) / "pipeline_info"
        pipeline_info_filename = next(pipeline_dir.glob("*_software_mqc_versions.yml"))

        with open(pipeline_info_filename) as f:
            return {
                software: version
                for software_dict in yaml.safe_load(f).values()
                for software, version in software_dict.items()
            }

    def get_run_parameters(self):
        results = OrderedDict()
        for key, value in self.run_parameters_tags.items():
            info = list(self.find(self.run_parameters, key))
            if info:
                results[value] = info[0]
        return results

    def get_read_cycles(self):
        read_and_cycles = OrderedDict()
        read_counter = 1
        index_counter = 1
        try:
            reads = self.run_info["RunInfo"]["Run"]["Reads"]["Read"]

            # Handle potential cases with only one Read element
            if isinstance(reads, dict):
                reads = [reads]  # Wrap it in a list to make it iterable

            for read_info in sorted(reads, key=itemgetter("@Number")):
                if read_info["@IsIndexedRead"] == "Y":
                    read_and_cycles[f"Index {index_counter} (bp)"] = read_info[
                        "@NumCycles"
                    ]
                    index_counter += 1
                else:
                    read_and_cycles[f"Read {read_counter} (bp)"] = read_info[
                        "@NumCycles"
                    ]
                    read_counter += 1
            return read_and_cycles

        except TypeError:
            return read_and_cycles

    def get_info(self):
        results = self.get_read_cycles()
        results.update(self.get_run_parameters())
        flowcell_type = self.find_flowcell_type_novaseqx()
        if flowcell_type:
            results.update(flowcell_type)

        return results

    def get_demultiplexing_info(self):
        try:
            return {"Demultiplexing": self.get_software_version(self.runfolder)}
        except FileNotFoundError:
            pass

        return {}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dumps a metadata yaml for MultiQC")
    parser.add_argument(
        "--runfolder", type=str, required=True, help="Path to runfolder"
    )

    args = parser.parse_args()
    runfolder = args.runfolder

    runfolder_info = RunfolderInfo(runfolder)
    info = runfolder_info.get_info()

    print(
        """
id: 'sequencing_metadata'
section_name: 'Sequencing Metadata'
plot_type: 'html'
description: 'regarding the sequencing run'
data: |
    <dl class="dl-horizontal">
"""
    )
    for k, v in info.items():
        print(f"        <dt>{k}</dt><dd><samp>{v}</samp></dd>")

    print("        <dt>Demultiplexing</dt>")
    for software, version in runfolder_info.get_demultiplexing_info().items():
        print(f"            <dd>{software}: <samp>{version}</samp></dd>")

    print("    </dl>")
