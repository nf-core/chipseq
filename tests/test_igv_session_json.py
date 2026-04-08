#!/usr/bin/env python3
"""
Tests for the IGV.js JSON session generation added to igv_files_to_session.py.

Run:  python -m pytest tests/test_igv_session_json.py -v
"""

import json
import os
import sys
import tempfile

# We can't import the script directly (it runs at module level) so we test via
# subprocess invocation — which is exactly how Nextflow calls it.
import subprocess

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "bin", "igv_files_to_session.py")


def _run(tmp_path, extra_args=None):
    """
    Create minimal input files, invoke the script, return the resulting
    session dict and the raw command output.
    """
    # -- list file: bigWig + narrowPeak entries (tab-separated: path\tcolour)
    list_file = os.path.join(tmp_path, "igv_files_orig.txt")
    with open(list_file, "w") as fh:
        fh.write("../../bwa/merged_library/bigwig/SPT5_T0_REP1.mLb.clN.bigWig\t0,0,178\n")
        fh.write("../../bwa/merged_library/bigwig/SPT5_T0_REP2.mLb.clN.bigWig\t0,0,178\n")
        fh.write("../../bwa/merged_library/macs3/narrow_peak/SPT5_T0_REP1_peaks.narrowPeak\t0,0,178\n")
        fh.write("../../bwa/merged_library/macs3/narrow_peak/SPT5_T0_REP2_peaks.narrowPeak\t0,0,178\n")

    replace_file = os.path.join(tmp_path, "replace_paths.txt")
    open(replace_file, "w").close()  # empty

    xml_out = os.path.join(tmp_path, "igv_session.xml")
    json_out = os.path.join(tmp_path, "igv_session.json")

    cmd = [
        sys.executable,
        SCRIPT,
        xml_out,
        list_file,
        replace_file,
        "../../genome/genome.fa",
        "--path_prefix", "../../",
        "--json_out", json_out,
        "--genome_id", "GRCh38",
        "--bam_files",
        "SPT5_T0_REP1.mLb.clN.sorted.bam",
        "SPT5_T0_REP2.mLb.clN.sorted.bam",
        "--bai_files",
        "SPT5_T0_REP1.mLb.clN.sorted.bam.bai",
        "SPT5_T0_REP2.mLb.clN.sorted.bam.bai",
        "--control_bam_files",
        "SPT5_INPUT.mLb.clN.sorted.bam",
        "--control_bai_files",
        "SPT5_INPUT.mLb.clN.sorted.bam.bai",
        "--sample_ids", "SPT5_T0_REP1", "SPT5_T0_REP2",
        "--control_ids", "SPT5_INPUT",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, f"Script failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"

    with open(json_out) as fh:
        session = json.load(fh)

    return session, result


def test_genome_is_resolved():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    assert session["genome"] == "hg38"


def test_genome_always_present():
    """genome key must always be present — Data Explorer requires it."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    assert "genome" in session


def test_genome_defaults_to_hg38_when_empty():
    """When --genome_id is empty, genome should default to hg38."""
    with tempfile.TemporaryDirectory() as tmp:
        list_file = os.path.join(tmp, "igv_files_orig.txt")
        with open(list_file, "w") as fh:
            fh.write("../../bwa/merged_library/bigwig/X.mLb.clN.bigWig\t0,0,178\n")
        replace_file = os.path.join(tmp, "replace_paths.txt")
        open(replace_file, "w").close()
        xml_out = os.path.join(tmp, "igv_session.xml")
        json_out = os.path.join(tmp, "igv_session.json")

        cmd = [
            sys.executable, SCRIPT, xml_out, list_file, replace_file,
            "../../genome/genome.fa",
            "--path_prefix", "../../",
            "--json_out", json_out,
            "--genome_id", "",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp)
        assert result.returncode == 0

        with open(json_out) as fh:
            session = json.load(fh)
        assert session["genome"] == "hg38"


def test_correct_number_of_tracks():
    """2 samples × 3 tracks + 1 control = 7 tracks total."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    assert len(session["tracks"]) == 7


def test_tracks_have_only_name_and_url():
    """Data Explorer schema: tracks only need name + url."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    for t in session["tracks"]:
        assert "name" in t
        assert "url" in t
        # Should NOT have type, format, color, height, order, indexURL etc.
        for forbidden in ("type", "format", "color", "height", "order", "indexURL", "autoscale", "displayMode"):
            assert forbidden not in t, f"track should not have '{forbidden}': {t}"


def test_urls_are_relative_no_dotdot():
    """JSON URLs should be relative from the output root (no ../../)."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    for t in session["tracks"]:
        assert not t["url"].startswith("../"), f"url still has ../ prefix: {t['url']}"


def test_control_track_present():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    ctrl = [t for t in session["tracks"] if "Input Control" in t["name"]]
    assert len(ctrl) == 1


def test_sample_names_in_track_names():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    names = [t["name"] for t in session["tracks"]]
    assert any("SPT5_T0_REP1" in n for n in names)
    assert any("SPT5_T0_REP2" in n for n in names)
    assert any("SPT5_INPUT" in n for n in names)


def test_xml_still_produced():
    """The original XML output should still be generated."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
        xml_path = os.path.join(tmp, "igv_session.xml")
        assert os.path.exists(xml_path)
        with open(xml_path) as fh:
            content = fh.read()
        assert '<?xml version="1.0"' in content
        assert "<Session" in content


def test_no_json_when_flag_omitted():
    """When --json_out is not provided, no JSON should be produced."""
    with tempfile.TemporaryDirectory() as tmp:
        list_file = os.path.join(tmp, "igv_files_orig.txt")
        with open(list_file, "w") as fh:
            fh.write("../../bwa/merged_library/bigwig/X.mLb.clN.bigWig\t0,0,178\n")
        replace_file = os.path.join(tmp, "replace_paths.txt")
        open(replace_file, "w").close()
        xml_out = os.path.join(tmp, "igv_session.xml")

        cmd = [sys.executable, SCRIPT, xml_out, list_file, replace_file, "hg38"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp)
        assert result.returncode == 0
        assert not os.path.exists(os.path.join(tmp, "igv_session.json"))


def test_broad_peak_url():
    """broadPeak files should be included as tracks with correct URL."""
    with tempfile.TemporaryDirectory() as tmp:
        list_file = os.path.join(tmp, "igv_files_orig.txt")
        with open(list_file, "w") as fh:
            fh.write("../../bwa/merged_library/bigwig/S1.mLb.clN.bigWig\t0,0,178\n")
            fh.write("../../bwa/merged_library/macs3/broad_peak/S1_peaks.broadPeak\t0,0,178\n")
        replace_file = os.path.join(tmp, "replace_paths.txt")
        open(replace_file, "w").close()
        xml_out = os.path.join(tmp, "igv_session.xml")
        json_out = os.path.join(tmp, "igv_session.json")

        cmd = [
            sys.executable, SCRIPT, xml_out, list_file, replace_file,
            "../../genome/genome.fa",
            "--path_prefix", "../../",
            "--json_out", json_out,
            "--genome_id", "GRCm38",
            "--bam_files", "S1.mLb.clN.sorted.bam",
            "--bai_files", "S1.mLb.clN.sorted.bam.bai",
            "--sample_ids", "S1",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp)
        assert result.returncode == 0

        with open(json_out) as fh:
            session = json.load(fh)

        assert session["genome"] == "mm10"
        peaks = [t for t in session["tracks"] if "Peaks" in t["name"]]
        assert len(peaks) == 1
        assert peaks[0]["url"].endswith(".broadPeak")
