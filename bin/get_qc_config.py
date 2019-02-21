from checkQC.run_type_recognizer import RunTypeRecognizer
from checkQC.config import ConfigFactory
import argparse

def convert_to_multiqc_config(checkqc_config):
    for qc_criteria in checkqc_config:
        if qc_criteria['name'] == 'ClusterPFHandler':
            bcl2fastq_total = {'total': {}}
            if not qc_criteria['warning'] == 'unknown':
                bcl2fastq_total['total']['warn'] = [{'lt': qc_criteria['warning']}]
            if not qc_criteria['error'] == 'unknown':
                bcl2fastq_total['total']['error'] = [{'lt': qc_criteria['error']}]

    return {'table_cond_formatting_rules': bcl2fastq_total}


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

    print(multiqc_config)
