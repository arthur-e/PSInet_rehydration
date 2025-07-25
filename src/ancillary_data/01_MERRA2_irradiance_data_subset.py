'''
Indexes down-welling, short-wave irradiance from MERRA-2 at the creosote site,
converts it to PAR.
'''

import datetime
import glob
import os
import re
import warnings
import netCDF4
import h5py
import numpy as np
import xarray as xr
import pandas as pd
from typing import Sequence
from functools import partial
from tqdm import tqdm

OUTPUT_CSV = '/home/arthur.endsley/Workspace/NTSG/projects/Y2026_PSInet/data/202050725_Jessica_Creosote_site_met_data_MERRA_SWGDN.csv'
EXPECTED_NUM_FILES = 91 # Daily files from 2023-04-01 through 2022-06-31 (N days)
MERRA2_DIR = '/anx_lagr4/SMAP/L4C_drivers/MERRA-2'
MERRA2_DATASETS = [ # (dataset, filename template) pairs
    ('tavg1_2d_rad_Nx',  'MERRA2_???.tavg1_2d_rad_Nx.*.???.nc*'),
]
MERRA2_FIELDS = {
    'tavg1_2d_rad_Nx':  ('SWGDN',)
}
COORDS = (-112.105, 32.753) # NOTE: Approximate
START_DATE = datetime.datetime(2023, 4, 1)
END_DATE = datetime.datetime(2023, 6, 30)

# Flip data up/down (ud) or left/right?
#   Argument lats/ lons should be a 1D numpy.ndarray
flip_ud = lambda lats: np.argmax(lats) > np.argmin(lats)
flip_lr = lambda lons: np.argmax(lons) < np.argmin(lons)


def main():
    'Populates the calibration driver dataset with data from MERRA2'
    # Get a list of dates in YYYYMMDD format
    date_list = [START_DATE]
    while date_list[-1] < END_DATE:
        date_list.append(date_list[-1] + datetime.timedelta(days = 1))
    date_list = [d.strftime('%Y%m%d') for d in date_list]
    coord_arrays = (
        xr.DataArray([COORDS[0]], dims = 'points'),
        xr.DataArray([COORDS[1]], dims = 'points'),
    )
    results = []
    for dataset, template in MERRA2_DATASETS:
        print(f'Working on "{dataset}" dataset...')
        # All of the following is to address the fact that MERRA-2 files have
        #   version numbers (UGH) in the filename that VARY between dates;
        #   hence, the filename is unpredictable
        # Get an initial list of possible files...
        some_list = glob.glob(f"{MERRA2_DIR}/{dataset.strip('*')}/{template}")
        some_dates = [f.split('.')[-3] for f in some_list]
        # some_list and some_date are NOT necessarily in chronological order;
        #   use our date_list to get the filenames in chron. order
        file_list = []
        for date in date_list:
            d = some_dates.index(date)
            file_list.append(some_list[d])
        assert len(file_list) == EXPECTED_NUM_FILES,\
            f'Did not find {EXPECTED_NUM_FILES} files!'
        # Iterate over files, in chronological order; "i" indexes the day
        for i, filename in enumerate(tqdm(file_list)):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                data = extract(filename, coord_arrays, MERRA2_FIELDS[dataset])
            for field_name in MERRA2_FIELDS[dataset]:
                results.append(pd.DataFrame({
                    'field': [field_name] * 24,
                    'date': [date_list[i]] * 24,
                    'hour': list(range(0, 24)),
                    'value': data[field_name].ravel()
                }))
        df = pd.concat(results)
        df.to_csv(OUTPUT_CSV)
        import ipdb
        ipdb.set_trace()#FIXME


def extract(filename: str, coords: tuple, fields: Sequence) -> list:
    '''
    Extracts data from a MERRA-2 NetCDF4 granule, based on the provided field.

    Parameters
    ----------
    filename : str
    coords : tuple
        A tuple of two floats, (longitude, latitude)
    feilds : Sequence
        Sequence of the field names to extract

    Returns
    -------
    dict
    '''
    nc = xr.open_dataset(filename)
    assert nc['time'].size == 24, 'Unexpected time axis length!'
    # Nearest-neighbour lookup for each site
    x, y = coords
    subset = nc.sel(lon = x, lat = y, method = 'nearest')
    results = []
    for field_name in fields:
        arr = subset[field_name].values
        results.append((field_name, arr))
    return dict(results)


if __name__ == '__main__':
    main()
