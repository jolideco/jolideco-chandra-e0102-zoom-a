import logging

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from gammapy.data import EventList
from gammapy.maps import Map

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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


def read_roi_geom(filename):
    """Read ROI geometry from counts map"""
    log.info(f"Reading {filename}")
    return Map.read(filename).geom


def select_events(filename, roi_geom):
    """Select events in ROI"""
    log.info(f"Reading {filename}")
    events = read_event_list_chandra(filename)
    selection = roi_geom.contains(events.radec)
    events_roi = events.table[selection]
    return events_roi


def bin_events(events, geom):
    """Bin events into counts map"""
    counts = Map.from_geom(geom)
    counts.fill_events(EventList(events))
    return counts


if __name__ == "__main__":
    filename_counts_in = snakemake.input.filename_counts
    filename_events_in = snakemake.input.filename_events

    events = read_event_list_chandra(filename_events_in)
    roi_geom = read_roi_geom(filename_counts_in)
    events_roi = select_events(events, roi_geom)

    iter_files = zip(snakemake.output.filename_events, snakemake.output.filename_counts)

    for filename_events, filename_counts in iter_files:
        random_state = np.random.RandomState(snakemake.wildcards.obs_id)

        events_bootstrap = Table(
            random_state.choice(events_roi, size=len(events_roi), replace=True)
        )
        events_bootstrap.write(filename_events, overwrite=True)

        counts = bin_events(events_bootstrap, roi_geom)
        counts.write(filename_counts, overwrite=True)
