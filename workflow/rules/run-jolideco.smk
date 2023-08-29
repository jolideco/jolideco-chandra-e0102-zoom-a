rule run_jolideco:
    input:
        expand("results/{{config_name}}/{obs_id}/maps/{{config_name}}-{obs_id}-counts.fits", obs_id=config["chandra-data"]["obs_ids"]),
        expand("results/{{config_name}}/{obs_id}/maps/{{config_name}}-{obs_id}-{irf_label}-{psf_simulator}-psf.fits", obs_id=config["chandra-data"]["obs_ids"], psf_simulator=config["chandra-data"]["psf_simulator"], irf_label=list(config["chandra-data"]["irfs"])),
        expand("results/{{config_name}}/{obs_id}/maps/{{config_name}}-{obs_id}-exposure.fits", obs_id=config["chandra-data"]["obs_ids"]),
    log:
        notebook="results/{config_name}/jolideco/{config_name}-jolideco.ipynb"
    output:
        filename_jolideco_result="results/{config_name}/jolideco/{config_name}-result-jolideco.fits",
        filename_npred_stacked="results/{config_name}/jolideco/{config_name}-npred.fits",
    # notebook:
    #     "../notebooks/jolideco-deconvolution.ipynb"
    shell:
        "cp {input[0]} {output.filename_jolideco_result}"
        "cp {input[0]} {output.filename_npred_stacked}"