nextflow.enable.dsl=2

include { HIFIASM } from './modules/nf-core/hifiasm/main'
include { YAK_COUNT } from './modules/nf-core/yak/count/main'


process BAM_FASTQ {
    tag "$meta.id"
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:df80ed8d23d0a2c43181a2b3dd1b39f2d00fab5c-0'
    input:
    tuple val(meta), path(bams)
    output:
    tuple val(meta), path("${meta.id}.fastq.gz")
    script:
    //def commands = bams.collect { "samtools quickcheck -u '${it}' && test \"\$(samtools view -c -F 4 '${it}')\" -eq 0 && samtools fastq -F 0x900 '${it}'" }.join(' && ')
    def commands = bams.collect { "samtools fastq -F 0x900 '${it}'" }.join(' && ')
    """
    set -o pipefail
    { $commands; } | pigz -p 4 -c > ${meta.id}.fastq.gz
    """
    stub:
    """
    printf '@read1\\nACGT\\n+\\nIIII\\n' | gzip -c > ${meta.id}.fastq.gz
    """
}

// Recover primary reads, including unmapped reads, singletons and duplicates.
process CRAM_FASTQ {
    tag "$meta.id"
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:df80ed8d23d0a2c43181a2b3dd1b39f2d00fab5c-0'
    input:
    tuple val(meta), path(cram), path(reference, name: 'reference.fa'), path(reference_fai, name: 'reference.fa.fai')
    output:
    tuple val(meta), path("${meta.id}.se.sr.fastq.gz")
    script:
    """
    set -o pipefail
    ${reference_fai ? '' : 'samtools faidx reference.fa'}
    samtools view -u -T reference.fa '${cram}' | samtools fastq -n -F 0x900 - | pigz -p 4 -c > '${meta.id}.se.sr.fastq.gz'
    """
    stub:
    """
    printf '@stub\\nACGT\\n+\\nIIII\\n' | gzip -c > '${meta.id}.se.sr.fastq.gz'
    """
}

process PIGZ_HIFIASM_OUTPUT {
    tag "$meta.id"
    cpus 4
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd226561e12b32bcacdf4f5ff74577e76233adf52ae5cbc499a2cdfe0e27d82/data'
    publishDir { "${params.outdir}/${meta.id}/hifiasm_extra" }, mode: 'copy'
    input:
    tuple val(meta), path(source)
    output:
    tuple val(meta), path("${source.name}.gz")
    script:
    """
    pigz -p ${task.cpus} -c '${source}' > '${source.name}.gz'
    """
    stub:
    """
    printf 'stub\\n' | gzip -c > '${source.name}.gz'
    """
}

process GFA_FASTA {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data'
    publishDir { "${params.outdir}/${meta.id}/assembly" }, mode: 'copy'
    input:
    tuple val(meta), path(primary), path(hap1), path(hap2)
    output:
    tuple val(meta), path('primary.fa'), path('hap1.fa'), path('hap2.fa')
    script:
    """
    for item in 'primary:${primary}' 'hap1:${hap1}' 'hap2:${hap2}'; do
        name=\${item%%:*}; gfa=\${item#*:}
        awk -F '\\t' -v h="\$name" '\$1 == "S" {if (\$3 == "*") exit 2; print ">" h "_" \$2 "\\n" \$3}' "\$gfa" > "\$name.fa"
        test -s "\$name.fa" || { echo "Empty \$name assembly" >&2; exit 1; }
    done
    """
    stub:
    """
    for name in primary hap1 hap2; do printf ">\${name}_ctg1\\nACGTACGT\\n" > \$name.fa; done
    """
}

process TELO_STATS {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2dae00435c9dff9202076069bf0d252db26a3e8c708741f17b5fbd053886ff9a/data'
    publishDir { "${params.outdir}/${meta.id}/qc" }, mode: 'copy'
    input:
    tuple val(meta), path(primary), path(hap1), path(hap2)
    output:
    tuple val(meta), path('*.telo.bed'), path('*.size.tsv')
    script:
    if (!(params.telo_motif ==~ /(?i)[ACGT]{1,15}/)) error 'telo_motif must contain 1–15 A/C/G/T bases'
    """
    for name in primary hap1 hap2; do
        seqtk telo -m ${params.telo_motif} \$name.fa > \$name.telo.bed
        seqtk size \$name.fa > \$name.size.tsv
    done
    """
    stub:
    """
    for name in primary hap1 hap2; do touch \$name.telo.bed; printf '1\\t8\\n' > \$name.size.tsv; done
    """
}

// Based on the nf-core gfastats process; grouped here to produce one sample tuple.
process ASSEMBLY_STATS {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7150110bea918bd523e41aab58f3be2eda949bbe7f715500681b01965e5f667/data'
    publishDir { "${params.outdir}/${meta.id}/qc" }, mode: 'copy'
    input:
    tuple val(meta), path(primary), path(hap1), path(hap2)
    output:
    tuple val(meta), path('*.gfastats.txt')
    script:
    """
    for name in primary hap1 hap2; do gfastats \$name.fa --threads ${task.cpus} > \$name.gfastats.txt; done
    """
    stub:
    """
    for name in primary hap1 hap2; do printf '# contigs: 1\\nTotal contig length: 8\\nContig N50: 8\\n' > \$name.gfastats.txt; done
    """
}

process ASMGENE_REFERENCE {
    tag 'reference transcripts'
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data'
    cpus 8
    memory '32 GB'
    input:
    path reference
    path transcripts
    output:
    path 'reference.cdna.paf'
    script:
    """
    minimap2 -cx splice:hq -t ${task.cpus} '${reference}' '${transcripts}' > reference.cdna.paf
    test -s reference.cdna.paf || { echo 'No reference transcript alignments' >&2; exit 1; }
    """
    stub:
    """
    printf 'stub\\n' > reference.cdna.paf
    """
}

process ASMGENE {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data'
    cpus 8
    memory '32 GB'
    publishDir { "${params.outdir}/${meta.id}/asmgene" }, mode: 'copy'
    input:
    tuple val(meta), path(primary), path(hap1), path(hap2)
    path refpaf
    path transcripts
    output:
    tuple val(meta), path('*.asmgene.tsv')
    script:
    """
    for name in primary hap1 hap2; do
        minimap2 -cx splice:hq -t ${task.cpus} \$name.fa '${transcripts}' > \$name.cdna.paf
        test -s \$name.cdna.paf || { echo "No transcript alignments for \$name" >&2; exit 1; }
        paftools.js asmgene '${refpaf}' \$name.cdna.paf > \$name.asmgene.tsv
    done
    """
    stub:
    """
    for name in primary hap1 hap2; do printf 'H\\tMetric\\treference\\tassembly\\nX\\tdup_cnt\\t10\\t9\\n' > \$name.asmgene.tsv; done
    """
}

process YAK_QV {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc0349f7b2bdd3d35e237ea545cb2523a91791cfa1246f4ac18c71aed303f2ab/data'
    cpus 8
    memory '32 GB'
    publishDir { "${params.outdir}/${meta.id}/qc" }, mode: 'copy'
    input:
    tuple val(meta), path(primary), path(hap1), path(hap2), path(yakdb)
    output:
    tuple val(meta), path('*.yak.qv.txt')
    script:
    """
    for name in primary hap1 hap2; do yak qv -t ${task.cpus} '${yakdb}' \$name.fa > \$name.yak.qv.txt; done
    """
    stub:
    """
    for name in primary hap1 hap2; do printf 'QV\\t40\\t40\\n' > \$name.yak.qv.txt; done
    """
}

process ALIGN_DIPLOID {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data'
    cpus 16
    memory '64 GB'
    input:
    tuple val(meta), path(hap1), path(hap2), path(reads)
    output:
    tuple val(meta), path('diploid.fa'), path('diploid.fa.fai'), path('reads.bam'), path('reads.bam.bai')
    script:
    """
    set -o pipefail
    cat hap1.fa hap2.fa > diploid.fa
    samtools faidx diploid.fa
    minimap2 -ax map-ont --eqx --cs -k17 -I8g -t ${task.cpus} diploid.fa '${reads}' | samtools sort -@ ${task.cpus} -o reads.bam
    samtools index -@ ${task.cpus} reads.bam
    """
    stub:
    """
    touch diploid.fa diploid.fa.fai reads.bam reads.bam.bai
    """
}

process FLAGGER {
    tag "$meta.id"
    container { params.flagger_sif ? file(params.flagger_sif).toAbsolutePath().toString() : "${projectDir}/assets/flagger.stub.sif" }
    cpus 16
    memory '64 GB'
    publishDir { "${params.outdir}/${meta.id}/flagger" }, mode: 'copy'
    input:
    tuple val(meta), path(assembly), path(fai), path(bam), path(bai)
    path alpha
    output:
    tuple val(meta), path('flagger_results')
    script:
    """
    awk 'BEGIN {OFS="\\t"} {print \$1,0,\$2}' '${fai}' > whole_genome.bed
    printf '{"whole_genome":"%s/whole_genome.bed"}\\n' "\$PWD" > annotations.json
    bam2cov --bam '${bam}' --output coverage.cov.gz --annotationJson annotations.json --threads ${task.cpus} --baselineAnnotation whole_genome
    mkdir flagger_results
    hmm_flagger --input coverage.cov.gz --outputDir flagger_results --preset ont-r10 --alphaTsv '${alpha}' --windowLen 8000 --labelNames Err,Dup,Hap,Col --threads ${task.cpus}
    test -s flagger_results/prediction_summary_final.tsv
    test -s flagger_results/final_flagger_prediction.bed
    """
    stub:
    """
    mkdir flagger_results
    printf 'stub; no biological results\\n' > flagger_results/prediction_summary_final.tsv
    """
}

process NUCFLAG {
    tag "$meta.id"
    container 'https://depot.galaxyproject.org/singularity/nucflag:1.0.2--pyhdfd78af_0'
    cpus 8
    memory '32 GB'
    publishDir { "${params.outdir}/${meta.id}/nucflag" }, mode: 'copy'
    input:
    tuple val(meta), path(assembly), path(fai), path(bam), path(bai)
    output:
    tuple val(meta), path('nucflag.bed')
    script:
    """
    nucflag call -i '${bam}' -f '${assembly}' -x ont_r10 -t ${task.cpus} -p 1 -o nucflag.bed
    """
    stub:
    """
    printf '# stub NucFlag\\n' > nucflag.bed
    """
}

process REPORT {
    tag "$meta.id"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy'
    input:
    tuple val(meta), path(telos), path(sizes), path(yak), path(flagger), path(stats), path(nucflag), path(asmgene)
    output:
    path 'summary.html'
    path 'summary.tsv'
    script:
    def opt = "${yak ? '--yak ' + yak.join(' ') : ''} ${flagger ? '--flagger ' + flagger : ''} --stats ${stats.join(' ')} ${nucflag ? '--nucflag ' + nucflag : ''} ${asmgene ? '--asmgene ' + asmgene.join(' ') : ''}"
    """
    assembly_report.py --sample '${meta.id}' --motif '${params.telo_motif}' --telo ${telos.join(' ')} --size ${sizes.join(' ')} $opt
    """
    stub:
    """
    printf 'Stub report: no biological results\\n' > summary.html
    printf 'sample\\tstatus\\n${meta.id}\\tstub\\n' > summary.tsv
    """
}

workflow {
    def cramFastqIntermediates = Collections.synchronizedSet([] as Set)
    workflow.onComplete {
        if (workflow.success) {
            cramFastqIntermediates.each { fastq ->
                try {
                    java.nio.file.Files.deleteIfExists(fastq)
                } catch (Exception e) {
                    log.warn "Could not remove intermediate ${fastq}: ${e.message}"
                }
            }
        }
    }
    if (!params.input) error 'Provide --input samples.csv'
    if (!params.skip_flagger && !workflow.stubRun) {
        if (!params.flagger_sif) error 'Provide --flagger_sif /shared/path/flagger-v1.2.0.sif (see README), or use --skip_flagger'
        def sif = file(params.flagger_sif, checkIfExists: true)
        if (!sif.toString().endsWith('.sif')) error 'flagger_sif must be a local Apptainer SIF file'
    }
    if ((!params.asmgene_reference) != (!params.asmgene_transcripts)) error 'Provide both --asmgene_reference and --asmgene_transcripts'
    if (!(params.yak_bloom.toString() ==~ /\d+/) || params.yak_bloom as int < 1 || params.yak_bloom as int > 40) error 'yak_bloom must be an integer from 1 to 40'
    def sheet = file(params.input, checkIfExists: true)
    def base = sheet.parent
    ch_samples = Channel.fromPath(sheet).splitCsv(header:true).map { row ->
        def id = row.sample?.trim()
        if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Invalid sample ID: ${id}"
        if (!row.bam?.trim()) error "Missing BAM for ${id}"
        def r1 = row.short_read_1?.trim() ?: ''
        def r2 = row.short_read_2?.trim() ?: ''
        def cram = row.short_read_cram?.trim() ?: ''
        def ref = row.cram_reference?.trim() ?: ''
        def fai = row.cram_reference_fai?.trim() ?: ''
        if ((!r1) != (!r2)) error "Both short reads must be supplied for ${id}"
        if ((!cram) != (!ref)) error "Both short_read_cram and cram_reference must be supplied for ${id}"
        if (fai && !ref) error "cram_reference_fai requires a CRAM reference for ${id}"
        if (cram && r1) error "Choose short-read CRAM or paired FASTQs, not both, for ${id}"
        if (r1 && (!(r1.endsWith('.gz')) || !(r2.endsWith('.gz')))) error "Short reads must be gzip-compressed FASTQ for ${id}"
        if (ref && ref.toLowerCase().endsWith('.gz')) error "CRAM reference must be an uncompressed FASTA for ${id}"
        def resolve = { x -> file(java.nio.file.Paths.get(x).isAbsolute() ? x : base.resolve(x), checkIfExists:true) }
        def refPath = ref ? resolve.call(ref) : ''
        def adjacentFai = ref ? file("${refPath}.fai") : null
        def faiPath = fai ? resolve.call(fai) : (adjacentFai?.exists() ? adjacentFai : '')
        tuple(id, resolve.call(row.bam.trim()), r1 ? resolve.call(r1) : '', r2 ? resolve.call(r2) : '', cram ? resolve.call(cram) : '', refPath, faiPath)
    }.groupTuple().map { id, bams, r1s, r2s, crams, refs, fais ->
        if ([r1s,r2s,crams,refs,fais].any { it.toSet().size()!=1 }) error "Short-read paths must match across BAM rows for ${id}"
        if (bams.toSet().size()!=bams.size()) error "Duplicate BAM for ${id}"
        tuple([id:id], bams, r1s[0], r2s[0], crams[0], refs[0], fais[0])
    }
    ch_samples = ch_samples.ifEmpty { error 'Samplesheet has no data rows' }
    ch_samples.multiMap { meta, bams, r1, r2, cram, ref, fai ->
        bam: tuple(meta,bams)
        shortreads: tuple(meta,r1,r2,cram,ref,fai)
    }.set { inputs }
    BAM_FASTQ(inputs.bam)
    HIFIASM(BAM_FASTQ.out.map { meta, reads -> tuple(meta, reads, []) }, Channel.value(tuple([:], [], [])), Channel.value(tuple([:], [], [])), Channel.value(tuple([:], [])))
    if (params.save_hifiasm_extra) {
        PIGZ_HIFIASM_OUTPUT(HIFIASM.out.corrected_reads.mix(HIFIASM.out.read_overlaps))
    }
    ch_primary = BAM_FASTQ.out.join(HIFIASM.out.primary_contigs, failOnMismatch:true).map { meta, reads, primary -> tuple(meta, primary) }
    ch_gfas = ch_primary.join(HIFIASM.out.hap1_contigs, failOnMismatch:true).join(HIFIASM.out.hap2_contigs, failOnMismatch:true)
    GFA_FASTA(ch_gfas)
    TELO_STATS(GFA_FASTA.out)
    ASSEMBLY_STATS(GFA_FASTA.out)
    if (params.asmgene_reference) {
        ch_transcripts = Channel.value(file(params.asmgene_transcripts, checkIfExists:true))
        ASMGENE_REFERENCE(file(params.asmgene_reference, checkIfExists:true), ch_transcripts)
        ASMGENE(GFA_FASTA.out, ASMGENE_REFERENCE.out, ch_transcripts)
        ch_asmgene = ASMGENE.out
    } else {
        ch_asmgene = inputs.bam.map { meta, bams -> tuple(meta,[]) }
    }
    ch_cram = inputs.shortreads.filter { meta, r1, r2, cram, ref, fai -> cram }.map { meta, r1, r2, cram, ref, fai -> tuple(meta, cram, ref, fai ? [fai] : []) }
    CRAM_FASTQ(ch_cram)
    ch_fastq = inputs.shortreads.filter { meta, r1, r2, cram, ref, fai -> r1 }.map { meta, r1, r2, cram, ref, fai -> tuple(meta, r1, r2) }
    ch_cram_fastq = CRAM_FASTQ.out.map { meta, fastq ->
        cramFastqIntermediates.add(fastq)
        tuple(meta + [single_end:true], fastq)
    }
    ch_short = ch_fastq.map { meta, r1, r2 -> tuple(meta + [single_end:false], [r1,r2]) }.mix(ch_cram_fastq)
    YAK_COUNT(ch_short)
    ch_yakdb = YAK_COUNT.out.yak.map { meta, db -> tuple([id:meta.id], db) }
    YAK_QV(GFA_FASTA.out.join(ch_yakdb))
    ch_flag_inputs = GFA_FASTA.out.join(BAM_FASTQ.out).filter { meta, primary, h1, h2, reads -> !params.skip_flagger || !params.skip_nucflag }.map { meta, primary, h1, h2, reads -> tuple(meta,h1,h2,reads) }
    ALIGN_DIPLOID(ch_flag_inputs)
    FLAGGER(ALIGN_DIPLOID.out.filter { meta, assembly, fai, bam, bai -> !params.skip_flagger }, file("${projectDir}/assets/ont-r10-alpha.tsv"))
    NUCFLAG(ALIGN_DIPLOID.out.filter { meta, assembly, fai, bam, bai -> !params.skip_nucflag })
    ch_nucflag = NUCFLAG.out.mix(inputs.bam.filter { meta,bams -> params.skip_nucflag }.map { meta,bams -> tuple(meta,[]) })
    ch_yak = YAK_QV.out.mix(inputs.shortreads.filter { meta,r1,r2,cram,ref,fai -> !r1 && !cram }.map { meta,r1,r2,cram,ref,fai -> tuple(meta,[]) })
    ch_flag = FLAGGER.out.mix(inputs.bam.filter { meta,bams -> params.skip_flagger }.map { meta,bams -> tuple(meta,[]) })
    REPORT(TELO_STATS.out.join(ch_yak).join(ch_flag).join(ASSEMBLY_STATS.out).join(ch_nucflag).join(ch_asmgene))
}
