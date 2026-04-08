#!/usr/bin/env python3

#######################################################################
#######################################################################
## Created on July 4th 2018 to create IGV session file from file list
## Extended 2025 to also produce an IGV.js JSON session for Seqera
## Data Explorer
#######################################################################
#######################################################################

import os
import errno
import argparse
import json
import re

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = (
    "Create IGV session file from a list of files and associated colours - "
    '".bed", ".bw", ".bigwig", ".tdf", ".gtf" files currently supported. '
    "Optionally emits an IGV.js-compatible JSON session."
)
Epilog = """Example usage: python igv_files_to_session.py <XML_OUT> <LIST_FILE> <REPLACE_FILE> <GENOME>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument("XML_OUT", help="XML output file.")
argParser.add_argument(
    "LIST_FILE",
    help="Tab-delimited file containing two columns i.e. file_name\\tcolour. Header isnt required.",
)
argParser.add_argument(
    "REPLACE_FILE",
    help="Tab-delimited file containing two columns i.e. file_name\\treplacement_file_name. Header isnt required.",
)
argParser.add_argument(
    "GENOME",
    help="Full path to genome fasta file or shorthand for genome available in IGV e.g. hg19.",
)

## OPTIONAL PARAMETERS
argParser.add_argument(
    "-pp",
    "--path_prefix",
    type=str,
    dest="PATH_PREFIX",
    default="",
    help="Path prefix to be added at beginning of all files in input list file.",
)
argParser.add_argument(
    "--json_out",
    type=str,
    dest="JSON_OUT",
    default="",
    help="If provided, also write an IGV.js-compatible JSON session to this path.",
)
argParser.add_argument(
    "--genome_id",
    type=str,
    dest="GENOME_ID",
    default="",
    help="Genome identifier for IGV.js (e.g. hg38, mm10). "
    "Used in the JSON output.  Falls back to GENOME if not set.",
)
argParser.add_argument(
    "--bam_files",
    nargs="*",
    default=[],
    help="BAM files to include in the JSON session (IP samples).",
)
argParser.add_argument(
    "--bai_files",
    nargs="*",
    default=[],
    help="BAM index files (same order as --bam_files).",
)
argParser.add_argument(
    "--control_bam_files",
    nargs="*",
    default=[],
    help="Input/control BAM files to include in the JSON session.",
)
argParser.add_argument(
    "--control_bai_files",
    nargs="*",
    default=[],
    help="Input/control BAM index files (same order as --control_bam_files).",
)
argParser.add_argument(
    "--sample_ids",
    nargs="*",
    default=[],
    help="Sample names for IP BAMs (same order as --bam_files). "
    "Derived from filenames when omitted.",
)
argParser.add_argument(
    "--control_ids",
    nargs="*",
    default=[],
    help="Sample names for control BAMs (same order as --control_bam_files). "
    "Derived from filenames when omitted.",
)
args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################


def makedir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


# ---- IGV.js genome name mapping -----------------------------------------

GENOME_MAP = {
    # Human
    "GRCh38": "hg38", "GRCh37": "hg19", "hg38": "hg38", "hg19": "hg19",
    "hg18": "hg18", "hs1": "hs1",
    # Mouse
    "GRCm39": "mm39", "GRCm38": "mm10", "mm39": "mm39", "mm10": "mm10",
    "mm9": "mm9",
    # Rat
    "Rnor_6.0": "rn6", "rn6": "rn6", "rn7": "rn7",
    # Fly / worm / fish / yeast
    "BDGP6": "dm6", "dm6": "dm6", "dm3": "dm3",
    "WBcel235": "ce11", "ce11": "ce11",
    "GRCz11": "danRer11", "GRCz10": "danRer10",
    "danRer11": "danRer11", "danRer10": "danRer10",
    "R64-1-1": "sacCer3", "sacCer3": "sacCer3",
    # Plant
    "TAIR10": "tair10", "tair10": "tair10",
    # Dog / cow
    "CanFam3.1": "canFam3", "canFam3": "canFam3",
    "UMD3.1": "bosTau8", "bosTau9": "bosTau9", "bosTau8": "bosTau8",
}


def resolve_genome_id(genome_key):
    """Return the IGV.js-compatible genome id for a given key."""
    return GENOME_MAP.get(genome_key, genome_key)


# ---- Per-sample colour palette (12 distinct, accessible colours) ---------

COLOUR_PALETTE = [
    "rgb(31,120,180)",   # blue
    "rgb(227,26,28)",    # red
    "rgb(51,160,44)",    # green
    "rgb(255,127,0)",    # orange
    "rgb(106,61,154)",   # purple
    "rgb(177,89,40)",    # brown
    "rgb(166,206,227)",  # light blue
    "rgb(251,154,153)",  # pink
    "rgb(178,223,138)",  # light green
    "rgb(253,191,111)",  # light orange
    "rgb(202,178,214)",  # light purple
    "rgb(255,255,153)",  # yellow
]

CONTROL_COLOUR = "rgb(128,128,128)"


def sample_name_from_bam(bam_path):
    """
    Derive a human-friendly sample name from a BAM / bigWig filename.

    nf-core/chipseq convention:  <sample>.mLb.clN.sorted.bam
    Falls back to the filename stem when the pattern doesn't match.
    """
    base = os.path.basename(bam_path)
    name = re.sub(r"\.mLb\..*$", "", base)
    name = re.sub(r"\.mLB\..*$", "", name)
    name = re.sub(r"\.sorted\.bam$", "", name)
    name = re.sub(r"\.bam$", "", name)
    name = re.sub(r"\.bigWig$", "", name)
    name = re.sub(r"\.bw$", "", name)
    return name


def sample_name_from_peak(peak_path):
    """Derive sample name from a MACS3 peak filename."""
    base = os.path.basename(peak_path)
    name = re.sub(r"_peaks\.(narrow|broad)Peak$", "", base)
    return name


############################################
############################################
## MAIN FUNCTION — XML (original)
############################################
############################################


def igv_files_to_session(XMLOut, ListFile, ReplaceFile, Genome, PathPrefix=""):
    makedir(os.path.dirname(XMLOut))

    replaceFileDict = {}
    fin = open(ReplaceFile, "r")
    while True:
        line = fin.readline()
        if line:
            ofile, rfile = line.strip().split("\t")
            replaceFileDict[ofile] = rfile
        else:
            break
            fin.close()
    fileList = []
    fin = open(ListFile, "r")
    while True:
        line = fin.readline()
        if line:
            ifile, colour = line.strip().split("\t")
            if len(colour.strip()) == 0:
                colour = "0,0,178"
            for ofile, rfile in replaceFileDict.items():
                if ofile in ifile:
                    ifile = ifile.replace(ofile, rfile)
            fileList.append((PathPrefix.strip() + ifile, colour))
        else:
            break
            fin.close()
    fout = open("igv_files.txt", "w")
    for ifile, colour in fileList:
        fout.write(ifile + "\n")
    fout.close()

    ## ADD RESOURCES SECTION
    XMLStr = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    XMLStr += '<Session genome="%s" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">\n' % (Genome)
    XMLStr += "\t<Resources>\n"
    for ifile, colour in fileList:
        XMLStr += '\t\t<Resource path="%s"/>\n' % (ifile)
    XMLStr += "\t</Resources>\n"

    ## ADD PANEL SECTION
    XMLStr += '\t<Panel height="1160" name="DataPanel" width="1897">\n'
    for ifile, colour in fileList:
        extension = os.path.splitext(ifile)[1].lower()
        if extension in [".bed", ".broadpeak", ".narrowpeak"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="10" height="20" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )
        elif extension in [".bw", ".bigwig", ".tdf"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="30" '
            XMLStr += (
                'id="%s" name="%s" normalize="false" renderer="BAR_CHART" sortable="true" visible="true" windowFunction="mean">\n'
                % (ifile, os.path.basename(ifile))
            )
            XMLStr += '\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10" minimum="0.0" type="LINEAR"/>\n'
            XMLStr += "\t\t</Track>\n"
        elif extension in [".gtf"]:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )
        elif extension in [".bam"]:
            pass
        else:
            XMLStr += (
                '\t\t<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="%s" '
                % (colour)
            )
            XMLStr += 'displayMode="SQUISHED" featureVisibilityWindow="-1" fontSize="10" height="20" '
            XMLStr += (
                'id="%s" name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n'
                % (ifile, os.path.basename(ifile))
            )
    XMLStr += "\t</Panel>\n"
    XMLStr += "</Session>"
    XMLOut = open(XMLOut, "w")
    XMLOut.write(XMLStr)
    XMLOut.close()

    return fileList


############################################
############################################
## IGV.js JSON SESSION BUILDER
############################################
############################################


def build_igvjs_session(
    genome_id,
    file_list,
    bam_files,
    bai_files,
    control_bam_files,
    control_bai_files,
    sample_ids,
    control_ids,
):
    """
    Build an IGV.js ``createBrowser()``-compatible configuration dict.

    Parameters
    ----------
    genome_id : str
        IGV.js genome identifier (e.g. "hg38").
    file_list : list[tuple[str, str]]
        (relative_path, colour) pairs already resolved by the XML builder
        for bigWigs, peaks, and consensus BEDs.
    bam_files / bai_files : list[str]
        IP / treatment BAM and index *filenames* (basenames — relative
        paths are constructed from the file_list prefix convention).
    control_bam_files / control_bai_files : list[str]
        Input / control BAM and index filenames (optional).
    sample_ids / control_ids : list[str]
        Explicit sample names; derived from filenames when empty.
    """
    igv_genome = resolve_genome_id(genome_id) if genome_id else "hg38"
    tracks = []

    # ------------------------------------------------------------------
    # Seqera Data Explorer renders IGV.js sessions with a minimal schema:
    #   { "genome": "...", "tracks": [{ "name": "...", "url": "..." }, ...] }
    #
    # All URLs are relative paths from the output directory root.  Data
    # Explorer resolves them against the bucket/directory containing the
    # JSON file.  Type and format are inferred from the file extension.
    # ------------------------------------------------------------------

    def _strip_dotdot(p):
        """Remove leading ../../ prefixes — give path relative to outdir root."""
        return re.sub(r"^(\.\./)+", "", p)

    # Determine aligner dir from existing paths (first bigWig or peak path)
    aligner_prefix = ""
    for fpath, _ in file_list:
        clean = _strip_dotdot(fpath)
        # e.g. "bwa/merged_library/bigwig/X.bigWig"
        parts = clean.split("/")
        if len(parts) >= 2:
            aligner_prefix = parts[0] + "/merged_library"
            break

    # ---- Group bigWig / peak paths by sample name ----
    bw_tracks = {}
    peak_tracks = {}
    consensus_tracks = []

    for fpath, colour in file_list:
        clean = _strip_dotdot(fpath)
        base = os.path.basename(fpath)
        ext = os.path.splitext(base)[1].lower()

        if ext in (".bw", ".bigwig"):
            name = sample_name_from_bam(base)
            bw_tracks[name] = (clean, colour)
            bw_tracks[name] = clean
        elif ext in (".narrowpeak", ".broadpeak"):
            name = sample_name_from_peak(base)
            peak_tracks[name] = clean
        elif ext == ".bed":
            consensus_tracks.append((clean, base))

    # ---- Determine sample ordering ----
    if sample_ids:
        ip_names = list(sample_ids)
    elif bam_files:
        ip_names = [sample_name_from_bam(b) for b in bam_files]
    else:
        # Fall back to names found in bigWig tracks
        ip_names = sorted(bw_tracks.keys())

    if control_ids:
        ctrl_names = list(control_ids)
    elif control_bam_files:
        ctrl_names = [sample_name_from_bam(b) for b in control_bam_files]
    else:
        ctrl_names = []

    # ---- Per-sample tracks ----
    for idx, sample in enumerate(ip_names):
        # BAM alignment — indexURL is required for IGV.js to load BAMs
        if idx < len(bam_files) and idx < len(bai_files):
            bam_base = os.path.basename(bam_files[idx])
            bai_base = os.path.basename(bai_files[idx])
            bam_url = "{}/{}".format(aligner_prefix, bam_base) if aligner_prefix else bam_base
            bai_url = "{}/{}".format(aligner_prefix, bai_base) if aligner_prefix else bai_base
            tracks.append({
                "name": "{} - Alignments".format(sample),
                "url": bam_url,
                "indexURL": bai_url,
            })

        # BigWig signal
        if sample in bw_tracks:
            tracks.append({
                "name": "{} - Signal".format(sample),
                "url": bw_tracks[sample],
            })

        # Peaks
        if sample in peak_tracks:
            tracks.append({
                "name": "{} - Peaks".format(sample),
                "url": peak_tracks[sample],
            })

    # ---- Control / input BAM tracks ----
    for idx in range(len(control_bam_files)):
        if idx >= len(control_bai_files):
            break
        ctrl_name = ctrl_names[idx] if idx < len(ctrl_names) else sample_name_from_bam(control_bam_files[idx])
        cbam_base = os.path.basename(control_bam_files[idx])
        cbai_base = os.path.basename(control_bai_files[idx])
        cbam_url = "{}/{}".format(aligner_prefix, cbam_base) if aligner_prefix else cbam_base
        cbai_url = "{}/{}".format(aligner_prefix, cbai_base) if aligner_prefix else cbai_base
        tracks.append({
            "name": "{} - Input Control".format(ctrl_name),
            "url": cbam_url,
            "indexURL": cbai_url,
        })

    # ---- Consensus peak tracks ----
    for clean_path, basename in consensus_tracks:
        tracks.append({
            "name": "{} - Consensus".format(os.path.splitext(basename)[0]),
            "url": clean_path,
        })

    session = {
        "genome": igv_genome,
        "tracks": tracks,
    }
    return session


############################################
############################################
## RUN FUNCTION
############################################
############################################

file_list = igv_files_to_session(
    XMLOut=args.XML_OUT,
    ListFile=args.LIST_FILE,
    ReplaceFile=args.REPLACE_FILE,
    Genome=args.GENOME,
    PathPrefix=args.PATH_PREFIX,
)

# ---- Optionally emit IGV.js JSON session ----
if args.JSON_OUT:
    genome_key = args.GENOME_ID if args.GENOME_ID else ""
    session = build_igvjs_session(
        genome_id=genome_key,
        file_list=file_list if file_list else [],
        bam_files=args.bam_files,
        bai_files=args.bai_files,
        control_bam_files=args.control_bam_files,
        control_bai_files=args.control_bai_files,
        sample_ids=args.sample_ids,
        control_ids=args.control_ids,
    )
    makedir(os.path.dirname(args.JSON_OUT))
    with open(args.JSON_OUT, "w") as fh:
        json.dump(session, fh, indent=2)
    print(
        "Wrote IGV.js session to {} (genome={}, {} tracks)".format(
            args.JSON_OUT, session["genome"], len(session["tracks"])
        )
    )
