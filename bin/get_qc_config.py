#!/usr/bin/env python

from checkQC.run_type_recognizer import RunTypeRecognizer
from checkQC.config import ConfigFactory
import argparse
import yaml


def get_mapping_handlers():
    multiqc_mapping = {'ClusterPFHandler':'total', #bcl2fastq
                       'ErrorRateHandler': 'Error', #Interop
                       'Q30Handler': 'percent_Q30', #bcl2fastq
                       'ReadsPerSampleHandler': 'mqc-generalstats-bcl2fastq-total'} #bcl2fastq

    return multiqc_mapping

def get_qc_criteria(mapper,checkqc_config):#This is ugly, checkq_config is passed through convert_to_multiqc_config. Maybe it should be global in a class.
    for handler in checkqc_config:
        if handler['name'] == mapper:
            return(handler)

def get_qc_level(mapper):
    #This is where we should get lt or gt.
    warning_level = {'ClusterPFHandler':'lt',
                     'ErrorRateHandler': 'gt',
                     'Q30Handler': 'lt',
                     'ReadsPerSampleHandler': 'lt'}

    return warning_level[mapper]


def convert_to_multiqc_config(checkqc_config):
    mapping_handlers = get_mapping_handlers()
    #I dont't want to return the whole dictionary from get_mapping_handlers, but I have to iterate over it.
    #Can I do this differently? It would be prefferd to iterate over the checkqc_config, but I don't want to use all keys
    #so then I have to make the code handle exceptions instead. Or should it be a global variable in an __init__ together with a config-file?
    multiqc_config_format = {}
    for mapper in mapping_handlers:
        multiqc_name = mapping_handlers[mapper]
        qc_criteria = get_qc_criteria(mapper,checkqc_config)
        qc_level = get_qc_level(mapper)
        multiqc_config_value = {multiqc_name: {}}
        if not qc_criteria['warning'] == 'unknown':
            multiqc_config_value[multiqc_name]['warn'] = [{qc_level: qc_criteria['warning']}]
        if not qc_criteria['error'] == 'unknown':
            multiqc_config_value[multiqc_name]['error'] = [{qc_level: qc_criteria['error']}]

        multiqc_config_format[multiqc_name] = multiqc_config_value[multiqc_name]

    return {'table_cond_formatting_rules': multiqc_config_format}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts CheckQC tresholds to MultiQC conditional format')
    parser.add_argument('--runfolder', type=str, required=True, help='Path to runfolder')
    parser.add_argument('--config', type=str, help='Path to checkQC config')

    args = parser.parse_args()
    runfolder = args.runfolder
    config = args.config

    run_type_recognizer = RunTypeRecognizer(runfolder)
    config = ConfigFactory.from_config_path(config)

    instrument_and_reagent_version = run_type_recognizer.instrument_and_reagent_version()
    both_read_lengths = run_type_recognizer.read_length()
    read_length = int(both_read_lengths.split("-")[0])
    checkqc_config = config.get_handler_configs(instrument_and_reagent_version, read_length)
    multiqc_config = convert_to_multiqc_config(checkqc_config)

    with open('qc_thresholds.yaml', 'w') as outfile:
        yaml.dump(multiqc_config, outfile)
