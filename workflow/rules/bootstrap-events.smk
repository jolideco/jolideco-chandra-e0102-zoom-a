rule bootstrap_events:
    input:
        filename_counts="results/{config_name}/{obs_id}/maps/{config_name}-{obs_id}-counts.fits",
        filename_events="results/{config_name}/{obs_id}/events/{config_name}-{obs_id}-events.fits",
    log:
        "logs/{config_name}/{obs_id}/bootstrap_events.log",
    output:
       filenames_events = expand("results/{{config_name}}/{{obs_id}}/events-bootstrapped/{{config_name}}-{{obs_id}}-iter-{idx}-events.fits",  idx=range(config["bootstrap"]["n_iter"])), 
       filenames_counts = expand("results/{{config_name}}/{{obs_id}}/maps-bootstrapped/{{config_name}}-{{obs_id}}-iter-{idx}-counts.fits", idx=range(config["bootstrap"]["n_iter"])) 

    script:
        "../scripts/bootstrap_events.py"