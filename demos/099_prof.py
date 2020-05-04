#creating paths
import sys
sys.path = ['/home/jluberto/kungpao'] + sys.path

#other general imports
import os
import warnings

import multiprocessing as mp
import numpy as np

from astropy.utils.console import ProgressBar
from astropy.io import fits

from crusher.galaxy import GalaxyMap
from crusher.data import BeneMassAgeZMaps
from crusher import visual
from crusher import combine
from crusher import profile

import math
import h5py

# HDF5 files
tng_file_hr = '/data/groups/leauthaud/simulations/TNG300/stellarmaps/galaxies_progmaps_hr_tng300_099.hdf5'
tng_file_lr = '/data/groups/leauthaud/simulations/TNG300/stellarmaps/galaxies_progmaps_lr_tng300_099.hdf5'

tng_label_hr = 'tng300_z0.4_hres'
tng_label_lr = 'tng300_z0.4_lres'

tng_data_hr = BeneMassAgeZMaps(tng_file_hr, label=tng_label_hr)
tng_data_lr = BeneMassAgeZMaps(tng_file_lr, label=tng_label_lr)

# making try, except because throws error if fits already created
try:
    tng_galaxies = tng_data_lr.sum_table(save=True)
    tng_galaxies.write(os.path.join(tng_data_lr.dir, '{}_galaxies.fits'.format(tng_data_lr.label)))
except:
    pass

# now for the code to process:
# gathering index for centrals/satellites and star mass for each galaxy
is_primary_hr = tng_data_hr.get_data(tng_file_hr)['catgrp_is_primary']
is_primary_lr = tng_data_lr.get_data(tng_file_lr)['catgrp_is_primary']

star_mass_hr = tng_data_hr.get_data(tng_file_hr)['scalar_star_mass']
star_mass_lr = tng_data_lr.get_data(tng_file_lr)['scalar_star_mass']

subhalo_pos_hr = tng_data_hr.get_data(tng_file_hr)['catsh_SubhaloCM']
subhalo_pos_lr = tng_data_lr.get_data(tng_file_lr)['catsh_SubhaloCM']

subhalo_id_hr = tng_data_hr.get_data(tng_file_hr)['catsh_id']
subhalo_id_lr = tng_data_hr.get_data(tng_file_lr)['catsh_id']

# problem galaxies in high res and low res (and first four)
prob_hr = (0, 1, 2, 3, 79, 497, 707, 866, 900, 1063, 1109, 1329, 1476, 2121, 2614, 3017, 3135, 3268, 3480, 3501, 3918, 4072)
prob_lr = (0, 1, 2, 3, 6, 49, 73, 182, 199, 214, 228, 326, 373, 386, 426, 487, 623, 664, 688, 739, 807, 893, 943, 975, 1337, 1532, 1563, 1651, 1739, 1809, 2019, 2032, 2184, 2223, 2230, 2234, 2385, 2418, 2454, 2501, 2561, 2608, 2665, 2836, 2903, 2966, 2970, 2976, 3016, 3174, 3381, 3408, 3719, 3745, 3748, 3758, 3834, 3894, 3984, 4039)

prob_id_hr = [subhalo_id_hr[i] for i in prob_hr]
prob_id_lr = [subhalo_id_lr[i] for i in prob_lr]

# function to get data from individual galaxies
def ell_data(tng_data, i, proj='xy'):
    gal = GalaxyMap(tng_data, i, proj=proj, aper_force=None)

    gal.aper_summary()
    gal.ell_summary()

    return visual.prepare_show_ellipse(gal.info, gal.maps, gal.ell_sum)['ell_gal_3']['intens']

profile_hr, profile_lr = list(), list()
for i in range(tng_data_lr.n_gal):

	if i in prob_hr:
		print('Dealing with galaxy:', i + 1)

		profile_hr.append(ell_data(tng_data_hr, i))

	if i in prob_lr:
		print('Dealing with galaxy:', i + 1)

		profile_lr.append(ell_data(tng_data_lr, i))

profile_hr, profile_lr = np.array(profile_hr), np.array(profile_lr)

# save to fits
hdu_hr = fits.PrimaryHDU(profile_hr)
hdu_lr = fits.PrimaryHDU(profile_lr)

hdul_hr = fits.HDUList([hdu_hr])
hdul_hr.writeto('099hr.fits')


hdul_lr = fits.HDUList([hdu_lr])
hdul_lr.writeto('099lr.fits')
