#!/usr/bin/env python3
import logging
import os
import sys

log = logging.getLogger(__name__)

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

name_template = "-".join(job_properties["wildcards"].values())
name = job_properties["rule"] + "-" + name_template

ARGS = {
    "-v": "JOLIDECO_GMM_LIBRARY",
    "-S": "/bin/bash",
    "-pe": "mthread 8",
    "-q": "sThC.q",
    "-l": "mres=64G,h_data=8G,h_vmem=8G",
    "-j": "y",
    "-N": f"{name}",
    "-M": "axel.donath@cfa.harvard.edu",
}

args = " ".join(f"{key} {value}" for key, value in ARGS.items())

cmd = f"qsub {args} {jobscript}"
log.info(cmd)
os.system(cmd)
