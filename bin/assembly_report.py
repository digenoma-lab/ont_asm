#!/usr/bin/env python3
import argparse
import csv
import html
from collections import Counter
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--sample', required=True)
p.add_argument('--motif', required=True)
p.add_argument('--telo', nargs='+', required=True)
p.add_argument('--size', nargs='+', required=True)
p.add_argument('--stats', nargs='+', required=True)
p.add_argument('--yak', nargs='*', default=[])
p.add_argument('--asmgene', nargs='*', default=[])
p.add_argument('--flagger')
p.add_argument('--nucflag')
p.add_argument('--readstats')
a = p.parse_args()

def by_assembly(paths):
    return {Path(path).name.split('.')[0]: Path(path) for path in paths}

def gfastats(path):
    fields = {}
    for line in path.read_text().splitlines():
        if ':' in line:
            key, value = line.split(':', 1)
            fields[key.strip()] = value.strip().replace(',', '')
    return fields

def asmgene(path):
    fields = {}
    for line in path.read_text().splitlines():
        item = line.split()
        if len(item) >= 4 and item[0] == 'X':
            fields[item[1]] = (item[2], item[3])
    return fields

def mmc(fields):
    try:
        reference, retained = map(int, fields['dup_cnt'])
        return f'{100 * (1 - retained / reference):.2f}' if reference else 'NA'
    except (KeyError, ValueError):
        return 'NA'

def readstats(path):
    if not path:
        return {}
    with Path(path).open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'), None) or {}
    return {key: value.replace(',', '') for key, value in row.items()}

sizes = by_assembly(a.size)
telos = by_assembly(a.telo)
stats = by_assembly(a.stats)
yaks = by_assembly(a.yak)
genes = by_assembly(a.asmgene)
rows = []
for name, size in sorted(sizes.items()):
    count, bases = size.read_text().split()[:2]
    ends = {}
    for line in telos[name].read_text().splitlines():
        fields = line.split('\t')
        if len(fields) < 4:
            continue
        contig, start, end, length = fields[:4]
        ends.setdefault(contig, set()).update((['left'] if int(start) == 0 else []) + (['right'] if int(end) == int(length) else []))
    qv = 'NA'
    if name in yaks:
        for line in yaks[name].read_text().splitlines():
            if line.startswith('QV\t'):
                qv = line.split('\t')[-1]
    gf = gfastats(stats[name])
    gene = asmgene(genes[name]) if name in genes else {}
    rows.append([a.sample, name, count, bases, gf.get('Contig N50', 'NA'), gf.get('Contig L50', 'NA'), gf.get('Largest contig', 'NA'), gf.get('GC content %', 'NA'), sum(len(x) == 2 for x in ends.values()), qv, mmc(gene)])
header = ['sample', 'assembly', 'contigs', 'bases', 'contig_n50', 'contig_l50', 'largest_contig', 'gc_percent', 'candidate_t2t_contigs', 'yak_adjusted_qv', 'asmgene_missing_multicopy_percent']
with open('summary.tsv', 'w') as out:
    w = csv.writer(out, delimiter='\t', lineterminator='\n'); w.writerow(header); w.writerows(rows)
flag = 'Not run'
if a.flagger:
    path = Path(a.flagger) / 'prediction_summary_final.tsv'
    flag = path.read_text() if path.exists() else 'Summary unavailable; inspect Flagger output.'
nuc = 'Not run'
if a.nucflag:
    counts = Counter()
    for line in Path(a.nucflag).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 4:
            continue
        hap = fields[0].split('_', 1)[0]
        counts[(hap, fields[3])] += 1
    nuc = '\n'.join(f'{hap}\t{label}\t{n}' for (hap, label), n in sorted(counts.items())) or 'No intervals reported.'
    nuc = 'haplotype\tclassification\tinterval_count\n' + nuc
esc = lambda s: html.escape(str(s))
table = ''.join('<tr>' + ''.join('<td>' + esc(c) + '</td>' for c in row) + '</tr>' for row in rows)
read = readstats(a.readstats)
read_metrics = [
    ('reads', 'num_seqs'),
    ('bases', 'sum_len'),
    ('read_n50', 'N50'),
    ('read_length_q1', 'Q1'),
    ('read_length_median', 'Q2'),
    ('read_length_q3', 'Q3'),
    ('mean_read_length', 'avg_len'),
    ('longest_read', 'max_len'),
    ('mean_quality', 'AvgQual'),
    ('q20_percent', 'Q20(%)'),
    ('q30_percent', 'Q30(%)'),
    ('gc_percent', 'GC(%)'),
]
read_section = ('<table><tr><th>metric</th><th>value</th></tr>' +
                ''.join('<tr><td>' + esc(name) + '</td><td>' + esc(read.get(key, 'NA')) + '</td></tr>' for name, key in read_metrics) +
                '</table><p>Read N50 and length quartiles are computed by seqkit on primary ONT reads extracted from the input BAMs. '
                'Quality columns summarise base qualities; they are not an assembly QV.</p>') if read else '<p>Not available.</p>'
gene_sections = ''.join('<h3>' + esc(name) + '</h3><pre>' + esc(path.read_text()) + '</pre>' for name, path in sorted(genes.items())) or '<p>Not run: provide a reference genome and matching transcripts.</p>'
Path('summary.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Assembly summary</title><style>body{font:16px sans-serif;margin:3em;max-width:1400px}table{border-collapse:collapse}td,th{padding:.5em;border:1px solid #aaa}pre{overflow:auto}</style></head><body><h1>Assembly summary: ' + esc(a.sample) + '</h1><table><tr>' + ''.join('<th>' + esc(x) + '</th>' for x in header) + '</tr>' + table + '</table><p>Candidate T2T: seqtk detects motif ' + esc(a.motif) + ' at both ends; this does not establish chromosome completeness. Haplotypes may be only partially phased.</p><h2>Read stats</h2>' + read_section + '<h2>asmgene</h2><p>Missing multi-copy percentage is relative to multi-copy genes in the supplied reference; genuine variation and transcript choice can affect it.</p>' + gene_sections + '<h2>NucFlag: combined haplotypes</h2><p>Counts are BED intervals by classification, including non-error classes; intervals may overlap and are not independent confirmed errors.</p><pre>' + esc(nuc) + '</pre><h2>HMM-Flagger</h2><pre>' + esc(flag) + '</pre></body></html>')
