""" Used to combine high res and low res TNG maps """

import numpy as np
import h5py

from scipy.interpolate import interp1d


# use two maps of after ell_summary() in crusher.galaxy as variables and return a combined map
def combine(highres_map, lowres_map):
  # first, get pixel scales
  pix_high = highres_map.pix
  pix_low = lowres_map.pix
  
  # now we want to interpolate the first and second maps
  r = np.linspace(0, 600 ** 0.25, 700 ** 0.25)
  
  interp_high = interp1d(highres_map['ell_gal_3']['r_kpc'] ** 0.25,
                         np.log10(highres_map['ell_gal_3']['intens']),
                         bounds_error=False, fill_value=np.nan,
                         kind='slinear')
  
  interp_low = interp1d(lowres_map['ell_gal_3']['r_kpc'] ** 0.25,
                         np.log10(lowres_map['ell_gal_3']['intens']),
                         bounds_error=False, fill_value=np.nan,
                         kind='slinear')
  
  # finish interpolating
  interp_high_map = interp_high(r)
  interp_low_map = interp_low(r)
  
  # split at 80kpc (want highres to be inner part of map and lowres to be outer part of map)
  mask_high = (r <= 80 ** 0.25)
  mask_low = (r > 80 ** 0.25)
  
  interp_high_map = interp_high_map[mask_high]
  interp_low_map = interp_low_map[mask_low]
  
  merged_map = np.append(interp_high_map, interp_low_map)
  
  return merged_map
