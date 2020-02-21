#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/chipseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'BamTools': ['v_bamtools.txt', r"bamtools (\S+)"],
    'deepTools': ['v_deeptools.txt', r"plotFingerprint (\S+)"],
    'Picard': ['v_picard.txt', r"([\d\.]+)-SNAPSHOT"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'Pysam': ['v_pysam.txt', r"(\S+)"],
    'MACS2': ['v_macs2.txt', r"macs2 (\S+)"],
    'HOMER': ['v_homer.txt', r"(\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}

results = OrderedDict()
results['nf-core/chipseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['BamTools'] = '<span style="color:#999999;\">N/A</span>'
results['deepTools'] = '<span style="color:#999999;\">N/A</span>'
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['Pysam'] = '<span style="color:#999999;\">N/A</span>'
results['MACS2'] = '<span style="color:#999999;\">N/A</span>'
results['HOMER'] = False
results['featureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/chipseq Software Versions'
section_href: 'https://github.com/nf-core/chipseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
