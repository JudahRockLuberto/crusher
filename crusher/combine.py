""" Used to combine high res and low res TNG maps """

import numpy as np

import h5py


# use two BeneMassAgeZMaps from crusher.data as variables and return a combined map
def combine(first_map, second_map):
  # first, get pixel scales
  pix_first = first_map.pix
  pix_second = second_map.pix
  
  
