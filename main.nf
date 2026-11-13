#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { summarize_mzML } from './module/summarize_mzML'
include { extract_MS1_matrix } from './module/extract_MS1_matrix'
include { plot_MS1_landscape } from './module/plot_MS1_landscape'
include { summarize_search_idXML } from './module/summarize_search_idXML'
include { summarize_filtered_idXML } from './module/summarize_filtered_idXML'
include { summarize_PSMFeatures } from './module/summarize_PSMFeatures'
include { summarize_ProteinGroups } from './module/summarize_ProteinGroups'
include { summarize_ProteinQuantifier } from './module/summarize_ProteinQuantifier'
include { collect_QCMetrics } from './module/collect_QCMetrics'

def resolveInput(def csv, String value) {
    if (!value?.trim()) error "Empty input path in ${csv}"
    def result = file(value.startsWith('/') ? value : "${csv.parent}/${value}")
    if (!result.isFile()) error "Missing input file: ${result}"
    return result
}

workflow {
    if (!(params.quantification_mode in ['TMT', 'LFQ']))
        error 'quantification_mode must be TMT or LFQ'
    if (params.quantification_mode == 'TMT' && !params.channel_csv)
        error 'TMT QC requires channel_csv'

    mzmlRows = Channel.fromPath(params.mzml_csv)
        .splitCsv(header: true)
        .toList()
        .map { [rows: it] }

    searchRows = Channel.fromPath(params.search_idxml_csv)
        .splitCsv(header: true)
        .toList()
        .map { [rows: it] }

    filteredRows = Channel.fromPath(params.filtered_idxml_csv)
        .splitCsv(header: true)
        .toList()
        .map { [rows: it] }

    manifests = mzmlRows.combine(searchRows)
        .combine(filteredRows)
        .map { mzPayload, rawPayload, finalPayload ->
            def mzRows = mzPayload.rows
            def rawRows = rawPayload.rows
            def finalRows = finalPayload.rows

            if (!mzRows || !rawRows || !finalRows) error 'All three manifests must contain rows'

            def safe = { value -> value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/ }
            def mz = [:]

            mzRows.each { row ->
                def sample = row.sample_name?.trim()
                def fraction = row.fraction_id?.trim()
                if (!safe(sample) || !safe(fraction) || mz.containsKey([sample, fraction]))
                    error "Invalid or duplicate mzML sample/fraction: ${sample}/${fraction}"
                mz[[sample, fraction]] = resolveInput(file(params.mzml_csv), row.path?.trim())
            }

            def runNames = mz.keySet().collect { key -> "${key[0]}_${key[1]}" }

            if ((runNames as Set).size() != runNames.size())
                error 'Sample/fraction names create duplicate mzML output filenames'

            def search = [:]
            rawRows.each { row ->
                def sample = row.sample_name?.trim()
                def engine = row.search_engine?.trim()
                if (!safe(sample) || !(engine in ['Comet', 'Sage', 'MSFragger']) ||
                    search.containsKey([sample, engine]))
                    error "Invalid or duplicate search sample/engine: ${sample}/${engine}"
                search[[sample, engine]] = resolveInput(file(params.search_idxml_csv), row.idXML?.trim())
            }

            def filtered = [:]
            finalRows.each { row ->
                def sample = (params.quantification_mode == 'TMT' ? row.plex_id :
                            row.sample_name)?.trim()
                if (!safe(sample) || filtered.containsKey(sample))
                    error "Invalid or duplicate filtered sample: ${sample}"
                filtered[sample] = resolveInput(file(params.filtered_idxml_csv), row.idXML?.trim())
            }

            def mzSamples = mz.keySet().collect { it[0] } as Set
            if (mzSamples != (search.keySet().collect { it[0] } as Set) ||
                mzSamples != (filtered.keySet() as Set))
                error 'mzML, search, and filtered manifests must describe the same samples'

            def mzEntries = mz.keySet().sort().collect { key -> tuple(key[0], key[1], mz[key]) }
            def rawEntries = search.keySet().sort().collect { key -> tuple(key[0], key[1], search[key]) }
            def finalEntries = filtered.keySet().sort().collect { sample ->
                def engines = search.keySet().findAll { it[0] == sample }.collect { it[1] }
                tuple(sample, engines.size() == 1 ? engines[0] : 'consensus', filtered[sample])
            }
            [mz: mzEntries, raw: rawEntries, filtered: finalEntries]
        }

    mzEntries = manifests.flatMap { it.mz }

    summarize_mzML(mzEntries)
    extract_MS1_matrix(mzEntries.map { sample, fraction, mzml ->
        tuple("${sample}_${fraction}", mzml)
    })
    plot_MS1_landscape(extract_MS1_matrix.out.matrix)
    summarize_search_idXML(manifests.flatMap { it.raw })
    summarize_filtered_idXML(manifests.flatMap { it.filtered })
    if (params.psm_features_csv) {
        featureRows = Channel.fromPath(params.psm_features_csv)
            .splitCsv(header: true)
            .toList().map { [rows: it] }

        featureEntries = featureRows.combine(manifests)
            .flatMap { rowPayload, manifest ->
                def safe = { value -> value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/ }
                def knownSearches = manifest.raw.collect { [it[0], it[1]] } as Set
                def knownFractions = manifest.mz.collect { [it[0], it[1]] } as Set
                def seen = [] as Set
                rowPayload.rows.collect { row ->
                    def sample = row.sample_name?.trim()
                    def engine = row.search_engine?.trim()
                    def route = row.route?.trim()
                    def fraction = row.fraction_id?.trim() ?: ''
                    if (!safe(sample) || !safe(engine) || !([sample, engine] in knownSearches) ||
                        !(route in ['MSBooster', 'PSMFeatureExtractor']) ||
                        (route == 'MSBooster' && !([sample, fraction] in knownFractions)) ||
                        (fraction && !safe(fraction)) || !seen.add([sample, engine, route, fraction]))
                        error "Invalid or duplicate PSM feature manifest row: ${row}"
                    def input = resolveInput(file(params.psm_features_csv), row.path?.trim())
                    if ((route == 'MSBooster' && !input.name.endsWith('.pin')) ||
                        (route == 'PSMFeatureExtractor' && !input.name.endsWith('.idXML')))
                        error "Wrong PSM feature file type for ${route}: ${input}"
                    tuple(sample, engine, route, fraction, input)
                }
            }
        summarize_PSMFeatures(featureEntries)
        featureMetrics = summarize_PSMFeatures.out.metrics
    } else {
        featureMetrics = Channel.empty()
    }
    summarize_ProteinGroups(file(params.protein_groups))

    def channelInput = file(params.channel_csv ?: params.mzml_csv)
    def assayInput = file(params.assay_mapping_csv ?: params.mzml_csv)
    quantInputs = Channel.of(
        tuple('protein', file(params.protein_quant_csv), file(params.mzml_csv),
              channelInput, assayInput),
        tuple('peptide', file(params.peptide_quant_csv), file(params.mzml_csv),
              channelInput, assayInput))

    summarize_ProteinQuantifier(quantInputs)

    allMetrics = summarize_mzML.out.metrics
        .mix(summarize_search_idXML.out.metrics)
        .mix(summarize_filtered_idXML.out.metrics)
        .mix(featureMetrics)
        .mix(summarize_ProteinGroups.out.metrics)
        .mix(summarize_ProteinQuantifier.out.metrics)
        .collect()

    collect_QCMetrics(allMetrics)
}
