rule prepare_spectrum:
    input:
        "resources/e0102-spectrum-plucinsky-2017/e0102-spectrum-plucinsky-2017-chart.dat",
    output:
        filename_spectrum="results/{config_name}/spectral-fit/{irf_label}/{config_name}-{irf_label}-source-flux-chart.dat"
    log:
        "results/{config_name}/{irf_label}-prepare-spectrum.log",
    localrule: True
    shell:
        "cp {input} {output.filename_spectrum}"