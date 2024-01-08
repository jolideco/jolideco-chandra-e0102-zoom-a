#!/usr/bin/env python
# coding: utf-8

# # Jolideco Analysis of Chandra Data of E0102

# In[ ]:


# import os

# os.environ.setdefault('JOLIDECO_GMM_LIBRARY', '/Users/adonath/github/jolideco/jolideco-gmm-prior-library')


# In[ ]:


# In[ ]:
import warnings
from itertools import zip_longest
from pathlib import Path

import matplotlib as mpl
import numpy as np
from astropy.visualization import simple_norm
from astropy.wcs.wcs import FITSFixedWarning
from gammapy.maps import Map
from jolideco.core import MAPDeconvolver, MAPDeconvolverResult
from jolideco.models import (
    FluxComponents,
    NPredCalibration,
    NPredCalibrations,
    NPredModels,
    SpatialFluxComponent,
)
from jolideco.priors import GaussianMixtureModel, GMMPatchPrior
from jolideco.utils.norms import IdentityImageNorm
from jolideco.utils.numpy import split_datasets_validation
from matplotlib import pyplot as plt

warnings.filterwarnings("ignore", category=FITSFixedWarning)


# In[ ]:


# # for vscode dark theme
# plt.style.use('dark_background')
# mpl.rcParams['figure.facecolor'] = '#25292E'


# In[ ]:


if "snakemake" in globals():
    filenames = [Path(_) for _ in snakemake.input]
    filename_jolideco_result = Path(snakemake.output.filename_jolideco_result)
    filename_npred_stacked = Path(snakemake.output.filename_npred_stacked)
    idx_iter = snakemake.wildcards.idx
else:
    idx_iter = 0
    config_name = "e0102-zoom-a"
    PATH_BASE = Path(f"../../results/{config_name}/")
    filenames = list(
        Path(f"../../results/{config_name}/").glob(
            f"*/maps-bootstrapped/{config_name}-*-{idx_iter}-counts.fits"
        )
    )
    filename_jolideco_result = (
        PATH_BASE / "jolideco" / f"{config_name}-result-jolideco.fits"
    )
    filename_npred_stacked = PATH_BASE / "jolideco" / "{config_name}-npred.fits"

obs_id_ref = 8365


# In[ ]:


filenames_counts = [_ for _ in filenames if "counts" in _.name]


# In[ ]:


datasets = {}


def read_dataset(filename_counts, idx):
    """Read counts, exposure and psf maps."""
    path = filename_counts.parent.parent / "maps"
    filename_exposure = path / filename_counts.name.replace(
        f"iter-{idx}-counts", "exposure"
    )
    filename_psf = path / filename_counts.name.replace(
        f"iter-{idx}-counts", "e0102-zoom-a-marx-psf"
    )
    counts = Map.read(filename_counts)

    psf = Map.read(filename_psf)
    psf = psf.cutout(psf.geom.center_skydir, width="15.8 arcsec")
    psf.data /= psf.data.sum()
    return {
        "counts": counts,
        "exposure": Map.read(filename_exposure),
        "psf": psf,
        "background": Map.from_geom(counts.geom) + 1e-2,
    }


for filename in filenames_counts:
    obs_id = filename.parts[-3]
    datasets[f"obs-id-{obs_id}"] = read_dataset(filename, idx=idx_iter)


# ## Counts

# In[ ]:


stacked = Map.from_geom(datasets[f"obs-id-{obs_id_ref}"]["counts"].geom)

for name, dataset in datasets.items():
    stacked += dataset["counts"]


# In[ ]:


stacked.plot(cmap="viridis", add_cbar=True)


# ## PSF

# In[ ]:


wcs = datasets[f"obs-id-{obs_id_ref}"]["psf"].geom.wcs

fig, axes = plt.subplots(
    ncols=5, nrows=5, subplot_kw={"projection": wcs}, figsize=(12, 12)
)

for ax, (name, dataset) in zip(axes.flat, datasets.items()):
    psf = dataset["psf"]
    psf.plot(ax=ax, cmap="viridis", add_cbar=True, stretch="log", vmin=0)
    ax.set_title(f"{name}")
    ax.axis("off")


# In[ ]:


def to_jolideco_dataset(maps, dtype=np.float32):
    """Convert Gammapy maps to Jolideco dataset."""
    return {
        "counts": maps["counts"].data.astype(dtype),
        "background": maps["background"].data.astype(dtype),
        "psf": {"e0102": maps["psf"].data.astype(dtype)},
        "exposure": maps["exposure"].data.astype(dtype),
    }


# In[ ]:


datasets_jolideco = {name: to_jolideco_dataset(maps) for name, maps in datasets.items()}


# ## Run Jolideco

# In[ ]:


gmm = GaussianMixtureModel.from_registry("jwst-cas-a-v0.1")
gmm.meta.stride = 4
print(gmm)


# In[ ]:


gmm.plot_mean_images(ncols=16, figsize=(12, 8))


# In[ ]:


patch_prior = GMMPatchPrior(
    gmm=gmm, cycle_spin=True, cycle_spin_subpix=True, norm=IdentityImageNorm()
)


shape = datasets_jolideco[f"obs-id-{obs_id_ref}"]["counts"].shape
flux_init = np.random.normal(loc=3, scale=0.01, size=shape).astype(np.float32)

component = SpatialFluxComponent.from_numpy(
    flux=flux_init,
    prior=patch_prior,
    use_log_flux=True,
    upsampling_factor=2,
)


components = FluxComponents()
components["e0102"] = component

print(components)


# In[ ]:


calibrations = NPredCalibrations()

for name in datasets_jolideco:
    calibration = NPredCalibration(background_norm=1.0, frozen=False)
    calibrations[name] = calibration


calibrations[f"obs-id-{obs_id_ref}"].shift_xy.requires_grad = False

print(calibrations)


# In[ ]:


deconvolve = MAPDeconvolver(n_epochs=10, learning_rate=0.1, beta=1.0)
print(deconvolve)


# In[ ]:


datasets_train = split_datasets_validation(datasets_jolideco, n_validation=0)


# In[ ]:


deconvolver = MAPDeconvolver(n_epochs=250, learning_rate=0.1, beta=1.0)


# In[ ]:


result = deconvolver.run(
    components=components,
    calibrations=calibrations,
    **datasets_train,
)

result.write(filename_jolideco_result, overwrite=True)


# In[ ]:


# result = MAPDeconvolverResult.read(filename_jolideco_result)


# In[ ]:


plt.figure(figsize=(12, 8))
result.plot_trace_loss()
plt.legend(loc="upper center", ncols=4)


# ## Results

# In[ ]:


counts = np.sum([_["counts"] for _ in datasets_jolideco.values()], axis=0)

fig, axes = plt.subplots(ncols=2, subplot_kw={"projection": wcs}, figsize=(14, 6))

norm_asinh = simple_norm(
    counts,
    min_cut=0,
    max_cut=2.5,
    stretch="power",
    power=1.0,
)


im = axes[0].imshow(counts, origin="lower", interpolation="None")
axes[0].set_title("Counts")
plt.colorbar(im)

im = axes[1].imshow(
    result.components.flux_upsampled_total_numpy,
    origin="lower",
    norm=norm_asinh,
    interpolation="gaussian",
)
axes[1].set_title("Deconvolved")
plt.colorbar(im)


# In[ ]:


print(calibrations)


# ## Residuals

# In[ ]:


geom = datasets[f"obs-id-{obs_id_ref}"]["counts"].geom


# In[ ]:


npreds = {}

for name, dataset in datasets_jolideco.items():
    model = NPredModels.from_dataset_numpy(
        dataset=dataset,
        components=result.components,
    )

    fluxes = result.components.to_flux_tuple()
    npred = model.evaluate(fluxes=fluxes).detach().numpy()[0, 0]
    npreds[name] = Map.from_geom(data=npred, geom=geom)


npreds_calibrated = {}

for name, dataset in datasets_jolideco.items():
    model = NPredModels.from_dataset_numpy(
        dataset=dataset, components=result.components, calibration=calibrations[name]
    )

    fluxes = result.components.to_flux_tuple()
    npred = model.evaluate(fluxes=fluxes).detach().numpy()[0, 0]
    npreds_calibrated[name] = Map.from_geom(data=npred, geom=geom)


# In[ ]:


npred_stacked = Map.from_geom(geom=geom)

for npred in npreds_calibrated.values():
    npred_stacked.stack(npred)


npred_stacked.write(filename_npred_stacked, overwrite=True)


# In[ ]:


fig, axes = plt.subplots(
    ncols=5,
    nrows=5,
    subplot_kw={"projection": wcs},
    gridspec_kw={"wspace": 0.2},
    figsize=(16, 16),
)


for name, ax in zip_longest(sorted(datasets_jolideco), axes.flat):
    if name is None:
        ax.set_visible(False)
        continue

    dataset = datasets[name]
    counts = dataset["counts"].sum_over_axes(keepdims=False).smooth(5)
    npred = npreds[name].smooth(5)

    residual = (counts - npred) / np.sqrt(npred)

    residual.plot(ax=ax, vmin=-0.5, vmax=0.5, cmap="RdBu", add_cbar=True)
    ax.set_title(f"{name}")


# In[ ]:


fig, axes = plt.subplots(
    ncols=5,
    nrows=5,
    subplot_kw={"projection": wcs},
    gridspec_kw={"wspace": 0.2},
    figsize=(16, 16),
)


for name, ax in zip_longest(sorted(datasets_jolideco), axes.flat):
    if name is None:
        ax.set_visible(False)
        continue

    dataset = datasets[name]
    counts = dataset["counts"].sum_over_axes(keepdims=False).smooth(5)
    npred = npreds_calibrated[name].smooth(5)

    residual = (counts - npred) / np.sqrt(npred)

    residual.plot(ax=ax, vmin=-0.5, vmax=0.5, cmap="RdBu", add_cbar=True)
    ax.set_title(f"{name}")


# In[ ]:


print(result.calibrations)
