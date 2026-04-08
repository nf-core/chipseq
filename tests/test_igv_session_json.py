#!/usr/bin/env python3
"""
Tests for the IGV.js JSON session generation added to igv_files_to_session.py.

Run:  python -m pytest tests/test_igv_session_json.py -v
"""

import json
import os
import sys
import tempfile

# Ensure bin/ is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

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


def test_correct_number_of_tracks():
    """2 samples × 3 tracks + 1 control = 7 tracks total."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    assert len(session["tracks"]) == 7


def test_track_types():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    types = [t["type"] for t in session["tracks"]]
    assert types.count("alignment") == 3  # 2 IP + 1 control
    assert types.count("wig") == 2
    assert types.count("annotation") == 2


def test_alignment_tracks_have_index_url():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    for t in session["tracks"]:
        if t["type"] == "alignment":
            assert "indexURL" in t
            assert t["indexURL"].endswith(".bai")


def test_urls_are_relative_no_dotdot():
    """JSON URLs should be relative from the output root (no ../../)."""
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    for t in session["tracks"]:
        assert not t["url"].startswith("../"), f"url still has ../ prefix: {t['url']}"
        if "indexURL" in t:
            assert not t["indexURL"].startswith("../"), f"indexURL still has ../ prefix: {t['indexURL']}"


def test_control_track_is_grey():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    ctrl = [t for t in session["tracks"] if "Input Control" in t["name"]]
    assert len(ctrl) == 1
    assert ctrl[0]["color"] == "rgb(128,128,128)"


def test_heights():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    for t in session["tracks"]:
        if t["type"] == "alignment" and "Input Control" not in t["name"]:
            assert t["height"] == 200
        elif t["type"] == "wig":
            assert t["height"] == 100
        elif t["type"] == "annotation":
            assert t["height"] == 50
        elif "Input Control" in t["name"]:
            assert t["height"] == 150


def test_distinct_colours_per_sample():
    with tempfile.TemporaryDirectory() as tmp:
        session, _ = _run(tmp)
    ip_tracks = [t for t in session["tracks"] if "Input Control" not in t["name"]]
    # Tracks for sample 1 should share a colour, different from sample 2
    rep1 = [t["color"] for t in ip_tracks if "REP1" in t["name"]]
    rep2 = [t["color"] for t in ip_tracks if "REP2" in t["name"]]
    assert len(set(rep1)) == 1, "REP1 tracks should share the same colour"
    assert len(set(rep2)) == 1, "REP2 tracks should share the same colour"
    assert rep1[0] != rep2[0], "Different samples should have different colours"


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


def test_broad_peak_format():
    """When input files are broadPeak, the format should be set correctly."""
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
        peaks = [t for t in session["tracks"] if t["type"] == "annotation"]
        assert len(peaks) == 1
        assert peaks[0]["format"] == "broadPeak"
