#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'NGI-ChIPseq': ['v_ngi_chipseq.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'Picard': ['v_picard.txt', r"version (\S+)"],
    'deepTools': ['v_deeptools.txt', r"plotFingerprint (\S+)"],
    'NGSplot': ['v_ngsplot.txt', r"Version: (\S+)"],
    'MACS2': ['v_macs2.txt', r"macs2 (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['NGI-ChIPseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['deepTools'] = '<span style="color:#999999;\">N/A</span>'
results['NGSplot] = '<span style="color:#999999;\">N/A</span>'
results['MACS2'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'ngi-chipseq'
section_name: 'NGI-ChIPseq Software Versions'
section_href: 'https://github.com/SciLifeLab/NGI-ChIPseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
