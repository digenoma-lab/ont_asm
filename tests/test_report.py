#!/usr/bin/env python3
"""Check report parsing against representative gfastats/asmgene/NucFlag output."""
import csv
import subprocess
import tempfile
from pathlib import Path

script = Path(__file__).resolve().parents[1] / 'bin' / 'assembly_report.py'
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    files = {
        'hap1.size.tsv': '2\t1200\n',
        'hap1.telo.bed': 'hap1_ctg1\t0\t100\t1000\nhap1_ctg1\t900\t1000\t1000\n',
        'hap1.gfastats.txt': '# contigs: 2\nContig N50: 1000\nContig L50: 1\nLargest contig: 1000\nGC content %: 41.8\n',
        'hap1.asmgene.tsv': 'H\tMetric\tref\tassembly\nX\tdup_cnt\t10\t9\n',
        'nucflag.bed': '#chrom\tstart\tend\tname\nhap1_ctg1\t1\t5\tmisjoin\n',
    }
    for name, data in files.items():
        (root / name).write_text(data)
    subprocess.run(['python3', str(script), '--sample', 's', '--motif', 'CCCTAA', '--telo', 'hap1.telo.bed', '--size', 'hap1.size.tsv', '--stats', 'hap1.gfastats.txt', '--asmgene', 'hap1.asmgene.tsv', '--nucflag', 'nucflag.bed'], cwd=root, check=True)
    with (root / 'summary.tsv').open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['contig_n50'] == '1000'
    assert row['candidate_t2t_contigs'] == '1'
    assert row['asmgene_missing_multicopy_percent'] == '10.00'
    assert 'hap1\tmisjoin\t1' in (root / 'summary.html').read_text()
