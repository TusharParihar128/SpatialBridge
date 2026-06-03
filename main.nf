#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Processes ─────────────────────────────────────────────

process QC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path h5_file
    path spatial_dir

    output:
    path "filtered_qc.h5ad"
    path "qc_plot.png"

    script:
    """
    python3 ${projectDir}/modules/qc/qc.py \
        ${h5_file} \
        ${spatial_dir} \
        .
    """
}

process NORMALIZE {
    publishDir "${params.outdir}/normalize", mode: 'copy'

    input:
    path filtered_h5ad

    output:
    path "normalized.h5ad"

    script:
    """
    python3 ${projectDir}/modules/normalize/normalize.py \
        ${filtered_h5ad} \
        .
    """
}

process CLUSTER {
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path normalized_h5ad

    output:
    path "clustered.h5ad"
    path "clusters_umap.png"

    script:
    """
    python3 ${projectDir}/modules/cluster/cluster.py \
        ${normalized_h5ad} \
        .
    """
}

process DECONV_KAGGLE {
    publishDir "${params.outdir}/deconv", mode: 'copy'

    output:
    path "cell_proportions.csv"

    script:
    """
    echo "Checking Kaggle status..."
    while true; do
        FULL=\$(kaggle kernels status ${params.kaggle_kernel} 2>&1)
        echo "\$FULL"
        if echo "\$FULL" | grep -q "COMPLETE"; then
            break
        elif echo "\$FULL" | grep -q "ERROR"; then
            echo "Kaggle job failed!"
            exit 1
        fi
        sleep 60
    done
    echo "Fetching results..."
    kaggle kernels output ${params.kaggle_kernel} -p .
    """
}

process TME {
    publishDir "${params.outdir}/tme", mode: 'copy'

    input:
    path clustered_h5ad
    path cell_proportions

    output:
    path "tme_results.csv"
    path "tme_plots.png"
    path "moransI.csv"

    script:
    """
    python3 ${projectDir}/modules/tme/tme.py \
        ${clustered_h5ad} \
        ${cell_proportions} \
        .
    """
}

// ── Workflow ──────────────────────────────────────────────

workflow {
    h5_file     = file(params.h5_file)
    spatial_dir = file(params.spatial_dir)
    ref_dir     = file(params.ref_dir)

    // Step 1: QC
    qc_out = QC(h5_file, spatial_dir)

    // Step 2: Normalize
    norm_out = NORMALIZE(qc_out[0])

    // Step 3: Cluster
    cluster_out = CLUSTER(norm_out)

    // Step 4: Deconvolution via Kaggle
    deconv_out = DECONV_KAGGLE()

    // Step 5: TME
    TME(cluster_out[0], deconv_out)
}
