rule bootstrap_events:
    input:
        filename_counts="results/{config_name}/{obs_id}/maps/{config_name}-{obs_id}-counts.fits",
        filename_events="results/{config_name}/{obs_id}/events/{config_name}-{obs_id}-events.fits",
    log:
        "logs/bootstrap_events.log",
    output:
       filenames_events = expand("results/{config_name}/{obs_id}/events-bootstrapped/{config_name}-{obs_id}-iter-{idx}-events.fits", config_name=config["config_name"], obs_id=config["chandra-data"]["obs_ids"], idx=range(N_ITER_EVENTS) 
       filenames_counts = expand("results/{config_name}/{obs_id}/events-bootstrapped/{config_name}-{obs_id}-iter-{idx}-events.fits", config_name=config["config_name"], obs_id=config["chandra-data"]["obs_ids"], idx=range(N_ITER_EVENTS) 
    script:
        "scripts/bootstrap_events.py"