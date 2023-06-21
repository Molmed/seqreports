#!/usr/bin/env python
from __future__ import print_function
import xmltodict
from collections import OrderedDict
import re
import argparse
import os
import json


class RunfolderInfo:
    def __init__(self, runfolder, bcl2fastq_outdir):
        self.runfolder = runfolder
        self.run_parameters = self.read_run_parameters()
        self.stats_json = self.read_stats_json(bcl2fastq_outdir)
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
            "Mode": "Flowcell type",
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

    def read_stats_json(self, bcl2fastq_outdir):
        stats_json_path = os.path.join(
            self.runfolder, bcl2fastq_outdir, "Stats/Stats.json"
        )
        if os.path.exists(stats_json_path):
            with open(stats_json_path) as f:
                return json.load(f)
        else:
            return None

    def get_bcl2fastq_version(self, runfolder):
        with open(os.path.join(runfolder, "bcl2fastq_version")) as f:
            bcl2fastq_str = f.read()
        return bcl2fastq_str.split("v")[1].strip()

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
            for read_info in self.stats_json["ReadInfosForLanes"][0]["ReadInfos"]:
                if read_info["IsIndexedRead"]:
                    read_and_cycles[f"Index {index_counter} (bp)"] = read_info[
                        "NumCycles"
                    ]
                    index_counter += 1
                else:
                    read_and_cycles[f"Read {read_counter} (bp)"] = read_info[
                        "NumCycles"
                    ]
                    read_counter += 1
            return read_and_cycles
        except TypeError:
            return read_and_cycles

    def get_info(self):
        results = self.get_read_cycles()
        results.update(self.get_run_parameters())
        if os.path.exists(os.path.join(self.runfolder, "bcl2fastq_version")):
            results["bcl2fastq version"] = self.get_bcl2fastq_version(self.runfolder)
        return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dumps a metadata yaml for MultiQC")
    parser.add_argument(
        "--runfolder", type=str, required=True, help="Path to runfolder"
    )
    parser.add_argument(
        "--bcl2fastq-outdir",
        type=str,
        default="Data/Intensities/BaseCalls",
        help="Path to bcl2fastq output folder relative to the runfolder",
    )

    args = parser.parse_args()
    runfolder = args.runfolder
    bcl2fastq_outdir = args.bcl2fastq_outdir

    runfolder_info = RunfolderInfo(runfolder, bcl2fastq_outdir)
    results = runfolder_info.get_info()

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
    for k, v in results.items():
        print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
    print("    </dl>")
