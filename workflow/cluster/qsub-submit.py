#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

name_template = "-".join(job_properties["wildcards"].values())
name = job_properties["rule"] + "-" + name_template

ARGS = {
    "-S": "/bin/bash",
    "-pe": "mthread 16",
    "-q": "sThC.q",
    "-l": "mres=64G,h_data=4G,h_vmem=4G",
    "-j": "y",
    "-N": f"{name}",
    "-M": "axel.donath@cfa.harvard.edu",
}

args = " ".join(f"{key} {value}" for key, value in ARGS.items())

os.system(f"qsub {args} {jobscript}")
