#!/usr/bin/env python3

import os
import sys
import errno
import argparse
import re


def parse_args(args=None):
    Description = "Reformat nf-core/chipseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
    SPT5_T0,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,SPT5,1,SPT5_INPUT,1
    SPT5_T0,SRR1822154_1.fastq.gz,SRR1822154_2.fastq.gz,SPT5,2,SPT5_INPUT,2
    SPT5_INPUT,SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R1.fastq.gz,SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R2.fastq.gz,1,,,
    SPT5_INPUT,SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R1.fastq.gz,SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R2.fastq.gz,2,,,
    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/chipseq/samplesheet/v2.1/samplesheet_test.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 3
        HEADER = ["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        ## Check sample entries
        for line_number, line in enumerate(fin, start=1):
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                # Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        "Invalid number of columns (found = {}, minimum = {})!".format(len(lspl),len(HEADER)),
                        "Line {}".format(line_number),
                        line,
                    )
                num_cols = len([x for x in lspl[: len(HEADER)] if x])
                if num_cols < MIN_COLS:
                    print_error(
                        "Invalid number of populated columns (found = {}, minimum = {})!".format(num_cols,MIN_COLS),
                        "Line {}".format(line_number),
                        line,
                    )

                ## Check sample name entries
                sample, fastq_1, fastq_2, replicate, antibody, control, control_replicate = lspl[: len(HEADER)]
                if sample.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                    sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line {}".format(line_number), line)
                if not re.match(r"^[a-zA-Z0-9_.-]+$", sample):
                    print_error(
                        "Sample name contains invalid characters! Only alphanumeric characters, underscores, dots and dashes are allowed.",
                        "Line {}".format(line_number),
                        line,
                    )

                ## Check FastQ file extension
                for fastq in [fastq_1, fastq_2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line {}".format(line_number), line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line {}".format(line_number),
                                line,
                            )

                ## Check replicate column is integer
                if not replicate.isdecimal():
                    print_error("Replicate id not an integer!", "Line {}".format(line_number), line)
                    sys.exit(1)

                ## Check antibody and control columns have valid values
                if antibody:
                    if antibody.find(" ") != -1:
                        print(f"WARNING: Spaces have been replaced by underscores for antibody: {antibody}")
                        antibody = antibody.replace(" ", "_")
                    if not control:
                        print_error(
                            "Both antibody and control columns must be specified!",
                            "Line {}".format(line_number),
                            line,
                        )

                if control:
                    if control.find(" ") != -1:
                        print(f"WARNING: Spaces have been replaced by underscores for control: {control}")
                        control = control.replace(" ", "_")
                    if not control_replicate.isdecimal():
                        print_error("Control replicate id not an integer!", "Line {}".format(line_number), line)
                        sys.exit(1)
                    control = "{}_REP{}".format(control, control_replicate)
                    if not antibody:
                        print_error(
                            "Both antibody and control columns must be specified!",
                            "Line {}".format(line_number),
                            line,
                        )

                ## Auto-detect paired-end/single-end
                sample_info = []  ## [single_end, fastq_1, fastq_2, replicate, antibody, control]
                ## Paired-end short reads
                if sample and fastq_1 and fastq_2:
                    sample_info = ["0", fastq_1, fastq_2, replicate, antibody, control]
                ## Single-end short reads
                elif sample and fastq_1 and not fastq_2:
                    sample_info = ["1", fastq_1, fastq_2, replicate, antibody, control]
                else:
                    print_error("Invalid combination of columns provided!", "Line {}".format(line_number), line)

                ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, replicate, antibody, control ]]}
                replicate = int(replicate)
                sample_info = sample_info + lspl[len(HEADER) :]
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = {}
                if replicate not in sample_mapping_dict[sample]:
                    sample_mapping_dict[sample][replicate] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample][replicate]:
                        print_error("Samplesheet contains duplicate rows!", "Line {}".format(line_number), line)
                    else:
                        sample_mapping_dict[sample][replicate].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "single_end",
                        "fastq_1",
                        "fastq_2",
                        "replicate",
                        "antibody",
                        "control",
                    ]
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that replicate ids are in format 1..<num_replicates>
                uniq_rep_ids = sorted(list(set(sample_mapping_dict[sample].keys())))
                if len(uniq_rep_ids) != max(uniq_rep_ids) or 1 != min(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1..<num_replicates>!",
                        "Sample",
                        "{}, replicate ids: {}".format(sample, ",".join([str(x) for x in uniq_rep_ids])),
                    )
                    sys.exit(1)

                ## Check that multiple replicates are of the same datatype i.e. single-end / paired-end
                if not all(
                    x[0][0] == sample_mapping_dict[sample][1][0][0] for x in sample_mapping_dict[sample].values()
                ):
                    print_error(
                        f"Multiple replicates of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                for replicate in sorted(sample_mapping_dict[sample].keys()):
                    ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                    if not all(
                        x[0] == sample_mapping_dict[sample][replicate][0][0]
                        for x in sample_mapping_dict[sample][replicate]
                    ):
                        print_error(
                            f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                            "Sample",
                            sample,
                        )

                    for idx, val in enumerate(sample_mapping_dict[sample][replicate]):
                        control = "_REP".join(val[-1].split("_REP")[:-1])
                        control_replicate = val[-1].split("_REP")[-1]
                        if control and (
                            control not in sample_mapping_dict.keys()
                            or int(control_replicate) not in sample_mapping_dict[control].keys()
                        ):
                            print_error(
                                f"Control identifier and replicate has to match a provided sample identifier and replicate!",
                                "Control",
                                val[4],
                            )

                    ## Write to file
                    for idx in range(len(sample_mapping_dict[sample][replicate])):
                        fastq_files = sample_mapping_dict[sample][replicate][idx]
                        sample_id = "{}_REP{}_T{}".format(sample, replicate, idx + 1)
                        if len(fastq_files) == 1:
                            fout.write(",".join([sample_id] + fastq_files) + ",\n")
                        else:
                            fout.write(",".join([sample_id] + fastq_files) + "\n")

    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
