from pathlib import Path

from sherpa_contrib.chart import save_chart_spectrum

with Path("sherpa-model-e0102.py").open("r") as fh:
    exec(fh.read())


e_min = 0.5
e_max = 7
e_width = 0.01

dataspace1d(e_min, e_max, e_width)

# set_analysis(None, "energy", "rate", factor=1)

save_chart_spectrum(
    "e0102-spectrum-plucinsky-2017-chart.dat",
    elow=e_min,
    ehigh=e_max,
    ewidth=e_width,
    clobber=True,
)
