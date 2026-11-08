# ONT de novo assembly

A Nextflow DSL2 workflow for basecalled ONT R10 simplex reads in **unaligned BAM** files. It extracts FASTQ, runs hifiasm 0.25.0 in `--ont` mode, converts the primary and two partially phased assemblies to FASTA, and evaluates them with gfastats, seqtk, NucFlag, optional Yak and asmgene, and HMM-Flagger. It runs Singularity-format containers with Apptainer by default. This is an initial implementation: **assembly and QC have only been validated with stubs; do not use it for production human analyses until a real-data integration run has passed**.

## Inputs

Create a CSV with this header (paths are absolute or relative to the samplesheet):

```csv
sample,bam,short_read_cram,cram_reference,cram_reference_fai,short_read_1,short_read_2
HG001,/data/HG001/run1.bam,/data/HG001/illumina.cram,/refs/GRCh38.fa,/refs/GRCh38.fa.fai,,
HG001,/data/HG001/run2.bam,/data/HG001/illumina.cram,/refs/GRCh38.fa,/refs/GRCh38.fa.fai,,
HG002,/data/HG002/reads.bam,,,,,
HG003,/data/HG003/reads.bam,,,,/data/HG003/R1.fastq.gz,/data/HG003/R2.fastq.gz
```

Each row is an ONT BAM, and rows with the same sample ID are combined into one assembly. For optional Yak QC, provide either a short-read CRAM plus its exact alignment reference as an **uncompressed FASTA**, or paired gzip-compressed FASTQ files. Repeat the same short-read inputs on each row for a sample. Leave all short-read fields empty if unavailable. Only one short-read CRAM per sample is currently supported. Paths are relative to the samplesheet or absolute.

Optionally supply the matching FASTA index in `cram_reference_fai`. If omitted, the workflow uses an existing index alongside the reference at `<cram_reference>.fai`; if none exists, it builds one in the task work directory. The index must correspond to the supplied FASTA, including its byte layout.

The CRAM is decoded against the supplied reference and streamed into a single `${sample}.se.sr.fastq.gz` file for Yak; read-name collation is unnecessary. ONT and CRAM FASTQ compression use four pigz workers, and each extraction task requests five CPUs to also accommodate Samtools. Primary reads are retained regardless of pairing, mapping status, duplicate marking or QC-fail flag; secondary and supplementary alignments are excluded. Hard-clipped bases and reads missing from the CRAM cannot be recovered. Lossy qualities, duplicate handling, selective read retention or alignment-based filtering in the supplied CRAM can affect the k-mer spectrum and QV; use a CRAM that retains the original sample's high-quality reads whenever possible. The reference must match the one used to encode the CRAM, including contig sequences; a CRAM index is not needed. The short-read FASTQ is an intermediate in the Nextflow work directory and is deleted after successful workflow completion; it is retained after failure. Deleting it means that a later `-resume` may need to repeat CRAM extraction.

ONT BAMs must be basecalled, unmapped, contain primary reads and base qualities, and have distinct file basenames within a sample. BAMs with mapped records are rejected. Samples must have different IDs, restricted to letters, digits, `_`, `-`, and `.`. Use paths without whitespace or shell metacharacters. Inputs must be accessible from compute nodes and inside containers. The workflow does not basecall, demultiplex, assess read depth or filter by quality.

## Running

Requires Nextflow >=25.10, supported Java, and Apptainer available on the login and compute nodes. Most images are native Singularity SIFs fetched from Seqera Wave or Galaxy and cached by Nextflow; use a persistent shared Apptainer cache on the cluster. Despite `/docker/` appearing in some Wave HTTPS URLs, those URLs identify native SIF blobs.

Flagger v1.2.0 requires a local SIF available to all nodes. If your site does not already provide one, create it once on a compatible machine with network access and sufficient disk space:

```bash
mkdir -p /shared/YOUR_PROJECT/containers
scripts/prepare-flagger-sif.sh /shared/YOUR_PROJECT/containers/flagger-v1.2.0.sif
```

This uses `apptainer pull` to convert the upstream Flagger OCI image into a SIF and verifies the two executables; it does not require a Docker installation or daemon. The import still needs access to the upstream image registry. Alternatively use an existing verified Flagger v1.2.0 SIF. The helper has not yet been tested on Kutral. With Flagger enabled, pass its shared SIF path:


```bash
nextflow run main.nf -profile local --input samples.csv --flagger_sif /shared/YOUR_PROJECT/containers/flagger-v1.2.0.sif --outdir results -resume
nextflow run main.nf -profile slurm --input samples.csv --flagger_sif /shared/YOUR_PROJECT/containers/flagger-v1.2.0.sif --outdir results -resume
```

The `slurm` profile sets only the executor. Before running on Kutral, provide a site config with the actual account, queue/partition, limits, shared work directory, Apptainer cache and bind settings; this repository does not assume cluster-specific values. The default hifiasm request is 32 CPUs, 256 GB and 72 hours, which is a starting point, not a verified human sizing. Configure processes in an extra `-c kutral.config`, for example:

```groovy
process {
  executor = 'slurm'
  queue = 'YOUR_PARTITION'
  clusterOptions = '--account=YOUR_ACCOUNT'
  withName: HIFIASM { cpus = 48; memory = '512 GB'; time = '96h' }
}
apptainer.cacheDir = '/shared/YOUR_PROJECT/apptainer-cache'
```

Run with `-c kutral.config -work-dir /shared/YOUR_PROJECT/work`. The default mapper uses `-I8g`, appropriate for a human-sized diploid reference; larger genomes may require increasing this value in `ALIGN_DIPLOID` to avoid a split index and unreliable mapping qualities. Mapping and sorting run simultaneously and may use more threads than their nominal allocation; tune resources before cluster use.

Useful options:

- `--skip_flagger`: skip the diploid coverage model and the Flagger SIF requirement, useful for haploid or bacterial tests.
- `--skip_nucflag`: skip NucFlag independently of Flagger. The read-to-diploid mapping is skipped only when both are disabled.
- `--asmgene_reference ref.fa --asmgene_transcripts transcripts.fa`: enable asmgene using a reference genome and its matching transcript FASTA. The reference transcript alignment is computed once and reused across samples. Use a reference and annotations appropriate for the organism; neither input is inferred from the CRAM reference. Without both inputs, asmgene is marked not run. Default asmgene thresholds apply, including sex chromosomes.
- `--hifiasm_args '-f0'`: add hifiasm arguments; `-f0` is recommended by hifiasm for small genomes. The `test` profile also supplies `-f0` and small resource requests. This workflow still requires all three hifiasm assembly graphs; some small/haploid datasets may not produce both haplotype graphs, and will fail clearly rather than report fabricated haplotypes.
- `--save_hifiasm_extra`: compress and publish corrected reads and overlap files if hifiasm produces them. This option does not itself enable their generation; use appropriate `--hifiasm_args` when needed. By default these files are not compressed or published.
- `--telo_motif CCCTAA`: telomere repeat motif (default human motif; 1–15 A/C/G/T bases). `seqtk telo` default detection thresholds apply.
- `--yak_bloom 37`: Yak Bloom filter exponent; reduce it for small genomes. Yak QV compares the sample's short-read k-mers to each assembly; it is not a parental phasing evaluation.

Do not use the human motif and diploid Flagger model as biological QC for *E. coli*. A bacterial real-data run has not yet been performed. A suitable test needs real R10 FASTQ with qualities converted to an unaligned BAM, enough coverage, and confirmation that hifiasm produces the expected outputs.

## Output and interpretation

Each sample has `graphs/` (hifiasm GFA), `assembly/` (primary, hap1 and hap2 FASTA with unique prefixed contig IDs), `qc/` (gfastats assembly summaries, seqtk telomere BED and size tables; Yak QV files when short reads exist), `nucflag/` (NucFlag BED on the combined haplotypes), `asmgene/` (per-assembly comparisons when enabled), `flagger/` (Flagger predictions and statistics), optionally `hifiasm_extra/` (pigz-compressed corrected reads and overlaps), and `summary.tsv` / `summary.html`. The HTML includes the raw asmgene and Flagger summaries and NucFlag interval counts by haplotype and class; `candidate_t2t_contigs` counts contigs with seqtk motif calls at both ends. A motif at both ends is **not proof of a complete chromosome**. Contig N50, L50, largest contig and GC content come from gfastats. asmgene reports gene-focused recovery and missing multicopy percentage relative to the chosen reference, not genome-wide completeness or proof of misassembly. NG50, reference-based structural variant and phase-switch validation are not performed.

With ONT reads alone, hifiasm's hap1/hap2 are **partially phased**, not guaranteed chromosome-wide parental haplotypes. Flagger and NucFlag share the read alignment to the combined haplotype assembly; NucFlag uses the ONT R10 preset. Its classified intervals are candidates, not confirmed errors, and interval counts can overlap. Flagger maps reads to their combined assembly and reports coverage-model classes using the ONT-R10 settings and alpha matrix. The workflow uses whole-genome annotation only; it does not project satellite annotations or run conservative-call filtering, so coverage bias can affect results. Yak QV requires suitable high-quality short reads from the assembled sample; parental reads and trio phasing are not currently supported.

## Validation and provenance

To reproduce a wiring-only stub run, use `tests/stub.sh` with Nextflow in `PATH`. Stubs create artificial placeholders and disable container execution; they are **not** a test of container availability, command compatibility, scientific outputs, or real sequencing data. On the development machine, the two-sample stub and the skip-Flagger variant completed with Nextflow 26.10.5 and Java 21. The development machine was ARM64, while cluster architecture and image compatibility have not been checked.

The hifiasm and Yak count processes, and the grouped gfastats process, are adapted from [nf-core hifiasm](https://github.com/nf-core/modules/tree/master/modules/nf-core/hifiasm) and [nf-core yak/count](https://github.com/nf-core/modules/tree/master/modules/nf-core/yak/count) and [nf-core gfastats](https://github.com/nf-core/modules/tree/master/modules/nf-core/gfastats), retrieved 2026-11-05; stubs/version evaluation were adjusted. Samtools/pigz uses the nf-core dragmap multi-tool SIF; seqtk, minimap2/samtools and the Python reporting environment use nf-core SIF references; local workflow-specific processes handle their wiring. The ONT alpha matrix and method follow [HMM-Flagger v1.2.0](https://github.com/mobinasri/flagger/tree/v1.2.0); [hifiasm's instructions](https://github.com/chhylp123/hifiasm) describe ONT FASTQ input, the `--ont` option and the meaning of its partially phased outputs. Containers are versioned but not all are digest-pinned, and should be mirrored and verified for a production deployment.

NucFlag uses the Galaxy-hosted NucFlag 1.0.2 SIF; its [documentation](https://github.com/logsdon-lab/NucFlag) describes the ONT R10 preset. asmgene follows [Heng Li’s method](https://lh3.github.io/2020/12/25/evaluating-assembly-quality-with-asmgene) and requires matching reference and transcript data. The workflow wiring has been exercised with stubs, the report parser with representative synthetic outputs, and the NucFlag SIF inspected for its command-line interface. The SIF is amd64 and cannot run on the development ARM64 host; full SIF execution and biological outputs must be validated on Kutral.
