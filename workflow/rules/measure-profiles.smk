rule measure_profiles:
    input:
        filename_events_bootstrap=expand("results/{{config_name}}/{obs_id}/events-bootstrapped/{{config_name}}-{obs_id}-iter-{{idx}}-events.fits", obs_id=config["chandra-data"]["obs_ids"]),
        filename_events=expand("results/{{config_name}}/{obs_id}/events/{{config_name}}-{obs_id}-events.fits", obs_id=config["chandra-data"]["obs_ids"]),
        filename_jolideco="results/{config_name}/jolideco/iter-{idx}/{config_name}-iter-{idx}-result-jolideco.fits",
        filename_ref="results/e0102-zoom-a/1308/maps/e0102-zoom-a-1308-counts.fits",
    log:
        "logs/{config_name}/iter-{idx}-measure-profiles.log",
    localrule: True
    output:
        filename_jolideco_profile = "results/{config_name}/profiles/iter-{idx}/{config_name}-iter-{idx}-profile-jolideco.fits",
        filename_counts_profile = "results/{config_name}/profiles/iter-{idx}/{config_name}-iter-{idx}-profile-counts.fits",
        filename_counts_stacked = "results/{config_name}/counts-stacked/{config_name}-iter-{idx}-stacked-counts.fits",
    script:
        "../scripts/measure_profile.py"