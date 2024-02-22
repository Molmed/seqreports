#!/usr/bin/env python

from checkQC.run_type_recognizer import RunTypeRecognizer
from checkQC.config import ConfigFactory
import argparse
import yaml


class ValueHandlerMapper(object):
    def __init__(self, handler_name, multiqc_mapping, compare_direction):
        self.handler_name = handler_name
        self.multiqc_mapping = multiqc_mapping
        self.compare_direction = compare_direction


class HandlerMapper:
    def __init__(self):
        self._mapper_list = [
            ValueHandlerMapper(
                handler_name="ClusterPFHandler",
                multiqc_mapping="total",
                compare_direction="lt",
            ),
            ValueHandlerMapper(
                handler_name="ErrorRateHandler",
                multiqc_mapping="Error",
                compare_direction="gt",
            ),
            ValueHandlerMapper(
                handler_name="Q30Handler",
                multiqc_mapping="percent_Q30",
                compare_direction="lt",
            ),
            ValueHandlerMapper(
                handler_name="ReadsPerSampleHandler",
                multiqc_mapping="mqc-generalstats-bcl2fastq-total",
                compare_direction="lt",
            ),
        ]

        self.mapping = self._convert_to_mappings(self._mapper_list)

    def _convert_to_mappings(self, mapper_list):
        mapper_dict = {}
        for mapper in mapper_list:
            mapper_dict[mapper.handler_name] = mapper
        return mapper_dict


def convert_to_multiqc_config(checkqc_config_dict):
    multiqc_config_format = {}
    handler_mapper = HandlerMapper()
    for mapper_name, mapper in handler_mapper.mapping.items():
        qc_criteria = checkqc_config_dict.get(mapper.handler_name)
        multiqc_config_value = {mapper.multiqc_mapping: {}}
        if not qc_criteria["warning"] == "unknown":
            multiqc_config_value[mapper.multiqc_mapping]["warn"] = [
                {mapper.compare_direction: qc_criteria["warning"]}
            ]
        if not qc_criteria["error"] == "unknown":
            multiqc_config_value[mapper.multiqc_mapping]["fail"] = [
                {mapper.compare_direction: qc_criteria["error"]}
            ]

        multiqc_config_format[mapper.multiqc_mapping] = multiqc_config_value[
            mapper.multiqc_mapping
        ]

    return {"table_cond_formatting_rules": multiqc_config_format}


def convert_to_dict(checkqc_config):
    checkqc_config_dict = {}
    for qc_handler in checkqc_config:
        checkqc_config_dict[qc_handler["name"]] = qc_handler

    return checkqc_config_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts CheckQC tresholds to MultiQC conditional format"
    )
    parser.add_argument(
        "--runfolder", type=str, required=True, help="Path to runfolder"
    )
    parser.add_argument("--config", type=str, help="Path to checkQC config")

    args = parser.parse_args()
    runfolder = args.runfolder
    config = args.config

    run_type_recognizer = RunTypeRecognizer(runfolder)
    config = ConfigFactory.from_config_path(config)

    instrument_and_reagent_version = (
        run_type_recognizer.instrument_and_reagent_version()
    )
    both_read_lengths = run_type_recognizer.read_length()
    read_length = int(both_read_lengths.split("-")[0])
    checkqc_config = config.get_handler_configs(
        instrument_and_reagent_version, read_length,
        use_closest_read_length=True
    )
    checkqc_config_dict = convert_to_dict(checkqc_config)
    multiqc_config = convert_to_multiqc_config(checkqc_config_dict)

    with open("qc_thresholds.yaml", "w") as outfile:
        yaml.dump(multiqc_config, outfile)
