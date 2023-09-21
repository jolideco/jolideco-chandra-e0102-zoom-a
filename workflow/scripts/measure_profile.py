from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from gammapy.data import EventList
from gammapy.estimators import ImageProfileEstimator
from gammapy.maps import Map, WcsGeom
from regions import RectangleSkyRegion
from scipy.ndimage import gaussian_filter

center = SkyCoord(16.0172, -72.0340, unit="deg")

region = RectangleSkyRegion(
    center=center,
    width=5 * u.arcsec,
    height=0.6 * u.arcsec,
    angle=45 * u.deg,
)


GEOM_PROFILE = WcsGeom.create(
    skydir=region.center,
    width=(region.height, region.width),
    binsz=0.02 * u.arcsec,
    frame="galactic",
)

GEOM_PROFILE.wcs.wcs.crota = [45, 45]

GEOM_REF = Map.read(snakemake.input.filename_ref).geom.upsample(2)


def wcs_from_header_chandra(header, x_col=11):
    """Create WCS from event file header

    Parameters
    ----------
    header : `~astropy.io.fits.Header`
        FITS header

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        WCS object
    """
    y_col = x_col + 1
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [header[f"TCRPX{x_col}"], header[f"TCRPX{y_col}"]]
    wcs.wcs.cdelt = [header[f"TCDLT{x_col}"], header[f"TCDLT{y_col}"]]
    wcs.wcs.crval = [header[f"TCRVL{x_col}"], header[f"TCRVL{y_col}"]]
    wcs.wcs.ctype = [header[f"TCTYP{x_col}"], header[f"TCTYP{y_col}"]]
    return wcs


def read_event_list_chandra(filename, hdu="EVENTS"):
    """Read event list"""
    hdu = fits.open(filename)[hdu]

    table = Table.read(hdu)

    wcs = wcs_from_header_chandra(header=hdu.header)

    for colname in table.colnames:
        table.rename_column(colname, colname.upper())

    ra, dec = wcs.wcs_pix2world(table["X"], table["Y"], 1)
    table["RA"] = ra * u.deg
    table["DEC"] = dec * u.deg

    mjd = table.meta["MJDREF"]
    mjd_int = np.floor(mjd).astype(np.int64)
    table.meta["MJDREFI"] = mjd_int
    table.meta["MJDREFF"] = mjd - mjd_int
    table.meta["TIMESYS"] = "tt"  # TODO: not sure tt is correct here...
    return EventList(table)


def read_map(filename, hdu="E0102"):
    """Read map from FITS file"""
    with fits.open(filename) as hdulist:
        hdu = hdulist[hdu]
        data = hdu.data

    return Map.from_geom(GEOM_REF, data=data)


def measure_profile(image):
    """Measure profile"""
    est = ImageProfileEstimator(axis="lat", method="mean")
    return est.run(image)


def get_counts_image(filenames):
    """Measure profile"""

    counts = Map.from_geom(GEOM_PROFILE)

    for filename in filenames:
        events = EventList.read(filename)
        counts.fill_events(events)

    return counts


if __name__ == "__main__":
    snakemake = locals()["snakemake"]
    flux = read_map(snakemake.input.filename_jolideco)
    flux = flux.interp_to_geom(GEOM_PROFILE)
    profile_jolideco = measure_profile(image=flux)
    profile_jolideco.table.write(
        snakemake.output.filename_jolideco_profile, overwrite=True
    )

    filenames_events = list(snakemake.input.filename_events_bootstrap)
    counts = get_counts_image(filenames_events)

    counts.write(snakemake.output.filename_counts_stacked, overwrite=True)
    profile_counts = measure_profile(image=counts)
    profile_counts.table.write(snakemake.output.filename_counts_profile, overwrite=True)
