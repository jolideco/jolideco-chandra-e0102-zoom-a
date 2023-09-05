rule run_jolideco:
    input:
        expand("results/{{config_name}}/{obs_id}/maps-bootstrapped/{{config_name}}-{obs_id}-iter-{{idx}}-counts.fits", obs_id=config["chandra-data"]["obs_ids"]),
    log:
        notebook="results/{config_name}/jolideco/iter-{idx}/{config_name}-iter-{idx}-jolideco.ipynb"
    output:
        filename_jolideco_result="results/{config_name}/jolideco/iter-{idx}/{config_name}-iter-{idx}-result-jolideco.fits",
        filename_npred_stacked="results/{config_name}/jolideco/iter-{idx}/{config_name}-iter-{idx}-npred.fits",
    notebook:
        "../notebooks/jolideco-deconvolution.ipynb"
