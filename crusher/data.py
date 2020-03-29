#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Deal with the Illustris or TNG maps."""

import os

import h5py
import numpy as np

from astropy.table import Table, Column


__all__ = [
    'BeneMassAgeZMaps',
    'Z_SUN',
    ]

# Solar metallicity
Z_SUN = 0.0134


class BeneMassAgeZMaps(object):
    """Stellar mass, age, metallicity maps provided by Benedikt Diemer.

    Parameters
    ----------
    hdf5_file : str
        Location of the HDF5 file of the map.
    label : str
        Label of this map file. e.g. "tng100_z0.4_high_resolution".
        Default: None

    """

    def __init__(self, hdf5_file, label=None):
        # read in the file itself
        self.hdf5_file = hdf5_file
        directory, file_name = os.path.split(hdf5_file)
        self.dir = directory
        self.hdf5_name = file_name

        # open and read in the file
        self.data = h5py.File(self.hdf5_file, 'r')
        
        # grab or differentiate useful information.
        # keys for available data
        self.keys = list(self.data.keys())
        
        # get the attributes in the config group
        self.config_maps = self.data['config']
        self.config_keys = [key for key in self.config_maps.attrs.keys()]
        
        #array of all values in file
        self.hdf5_values = list(self.data)

        # redshift of the snapshot
        self.redshift = self.get_config_key('snap_z')

        # pixel scale in kpc/pixel
        self.pix = self.get_pixel_scale()

        # label the dataset
        if label is not None:
            self.label = label
        else:
            self.label = "{}_z{}".format(self.get_config_key('sim'), self.redshift)

        # number of galaxies
        self.n_gal = len(self.data['catsh_id'])

    def get_pixel_scale(self):
        """Get the pixel scale of the map."""
        return 2. * self.get_config_key('map_range_min') / self.get_config_key('map_npixel')

    def get_config_key(self, key):
        """Get the key value from the map config."""
        if key in self.config_keys:
            return self.config_maps.attrs[key]
        else:
            print("# Key: {} is not available!".format(key))
            return None

    def get_data(self, hdf5_file):
        """Extract all useful data from the HDF5 file.

        Parameters
        ----------
        hdf5_file : str
            Location of the HDF5 file of the map.

        """
        
        # array of useful values in this array (note some TNG files have different names for the same thing: 
        # 'scalar_star_mass' and 'scalar_stelar_mass')
        useful_values = ['catsh_id', 'catgrp_is_primary', 'scalar_star_mass', 'catgrp_Group_M_Crit200',
                         'catsh_SubhaloCM', 'catsh_SubhaloPos', 'scalar_stellar_mass',
                         'scalar_star_age', 'scalar_star_metallicity', 'map_star_rho_insitu_xy', 
                         'map_star_rho_insitu_xz', 'map_star_rho_insitu_yz', 'map_star_rho_exsitu_xy',
                         'map_star_rho_exsitu_xz', 'map_star_rho_exsitu_yz', 'map_star_age_insitu_xy',
                         'map_star_age_insitu_xz', 'map_star_age_insitu_yz', 'map_star_age_exsitu_xy',
                         'map_star_age_exsitu_xz', 'map_star_age_exsitu_yz', 'map_star_metallicity_insitu_xy',
                         'map_star_metallicity_insitu_xz', 'map_star_metallicity_insitu_yz', 'map_star_metallicity_exsitu_xy',
                         'map_star_metallicity_exsitu_xz', 'map_star_metallicity_exsitu_yz']
        
        # add useful values to array
        data_dict = dict()
        for i in range(len(self.hdf5_values)):
            print(i)
            if i in useful_values:
                data_dict[i] = list(self.data[i])

        return data_dict

    def sum_table(self, save=False):
        """Put the basic information of all galaxies in a table.

        Parameters
        ----------
        save : bool, optional
            Save the summary table in a npy output. Default: False.

        """
        
        # array of wanted values in the table
        wanted_values = ['n_gal', 'catsh_id', 'catgrp_is_primary', 'scalar_star_mass', 'catgrp_Group_M_Crit200',
                         'scalar_star_age', 'scalar_star_metallicity']
        wanted_values_names = ['index', 'catsh_id', 'cen_flag', 'logms', 'logm200c', 'age', 'metallicity']
        
        #add into table
        summary = Table()
        for i in range(len(wanted_values)):
            summary.add_column(Column(data))
        
        
        # because individual columns require individual equations (log, divide by solar metalicity, etc),
        # we must check if in self.hdf5 values individually
                                     
        # add columns that -should- be in evert dataset                             
        summary.add_column(Column(data=np.arange(self.n_gal), name='index'))
        summary.add_column(Column(data=np.asarray(self.data['catsh_id']), name='catsh_id'))
                                      
        # add columns that -might- be in the data
        if 'catgrp_is_primary' in self.hdf5_values:
            summary.add_column(Column(data=np.asarray(self.data['catgrp_is_primary']), name='cen_flag'))
        # i've seen both 'scalar_star_mass' and 'scalar_stellar_mass,' so both are added (assuming
        # a file doesn't have both)
        if 'scalar_star_mass' in self.hdf5_values:
            summary.add_column(Column(data=np.log10(np.asarray(self.data['scalar_star_mass'])), name='logms'))
        if 'scalar_stellar_mass' in self.hdf5_values:
            summary.add_columns(Column(data=np.log10(np.asarray(self.data['scalar_stellar_mass'])), name='logms'))
        if 'catgrp_Group_M_Crit200' in self.hdf5_values:
            summary.add_column(Column(data=np.log10(np.asarray(self.data['catgrp_Group_M_Crit200'])), name='logm200c'))
        if 'scalar_star_age' in self.hdf5_values:
            summary.add_column(Column(data=np.asarray(self.data['scalar_star_age']), name='age'))
        if 'scalar_star_age' in self.hdf5_values:
            summary.add_column(Column(data=np.log10(np.asarray(self.data['scalar_star_metallicity']) / Z_SUN),
            name='metallicity'))

        if save:
            np.save(os.path.join(self.dir, "{}_galaxies.npy".format(self.label)), summary)

        return summary

    def get_basic_info(self, idx):
        """Gather basic information of the galaxy.

        Parameters
        ----------
        idx: int
            Index of the galaxy.

        Return
        ------
        info: dict
            A dictionary that contains basic information of the galaxy.
        """
        
        info_dict = dict()
        
        
        # add possible values in data
        if 'catgrp_is_primary' in self.hdf5_values:
            info_dict['cen_flag'] = self.data['catgrp_is_primary'][idx]
        # just like in sum_table and get_data, do both 'scalar_star_mass' and 'scalar_stellar_mass'
        if 'scalar_star_mass' in self.hdf5_values:
            info_dict['logms'] = np.log10(self.data['scalar_star_mass'][idx])
        if 'scalar_stellar_mass' in self.hdf5_values:
            info_dict['logms'] = np.log10(self.data['scalar_stellar_mass'][idx])
        if 'catgrp_Group_M_Crit200' in self.hdf5_values:
            info_dict['logm200c'] = np.log10(self.data['catgrp_Group_M_Crit200'][idx])
        if 'scalar_star_age' in self.hdf5_values:
            info_dict['age'] = self.data['scalar_star_age'][idx]
        if 'scalar_star_metallicity' in self.hdf5_values:
            info_dict['metallicity'] = self.data['scalar_star_metallicity'][idx]
            
        # add other values that -should- be in data
        info_dict['pix'] = self.pix
        info_dict['catsh_id'] = self.data['catsh_id'][idx]
        
        return info_dict

    def get_maps(self, idx, proj, verbose=False, maps_only=False):
        """Gather the stellar mass, age, metallicity map.

        Parameters
        ----------
        idx: int
            Index of the galaxy.
        proj: str
            Projection of the map. [xy|xz|yz]

        Return
        ------
        info: dict
            A dictionary that contains basic information of the galaxy.
        maps: dict
            A dictionary that contains all the necessary maps.

        """
        # maps dict to add to later
        maps = dict()
        
        
        # Basic information about the content
        info = self.get_basic_info(idx)
        if verbose:
            print("\n# Subhalo ID: {}".format(info['catsh_id']))

        # Projection
        if proj not in ['xy', 'xz', 'yz']:
            raise Exception("# Wrong projection: [xy | xz | yz]")
        info['proj'] = proj

        # Get the stellar mass maps if there
        if 'map_star_rho_insitu_{}'.format(proj) in self.hdf5_values and 'map_star_rho_exsitu_{}'.format(proj) in self.hdf5_values:
            mass_ins = self.data['map_star_rho_insitu_{}'.format(proj)][idx] * (self.pix ** 2)
            mass_exs = self.data['map_star_rho_exsitu_{}'.format(proj)][idx] * (self.pix ** 2)
            mass_gal = mass_ins + mass_exs
            
            # Stellar mass on the maps
            info['logms_map_ins'] = np.log10(mass_ins.sum())
            info['logms_map_exs'] = np.log10(mass_exs.sum())
            info['logms_map_gal'] = np.log10(mass_gal.sum())
            
            # Image size
            img_h, img_w = mass_ins.shape
            info['img_h'] = img_h
            info['img_w'] = img_w
            info['img_cen_x'] = img_h / 2.
            info['img_cen_y'] = img_w / 2.
            
            if verbose:
                print("\n# log(M*_ins): {:6.2f}".format(info['logms_map_ins']))
                print("# log(M*_exs): {:6.2f}".format(info['logms_map_exs']))
                print("# log(M*_gal): {:6.2f}".format(info['logms_map_gal']))
            
            maps['mass_ins'] = mass_ins
            maps['mass_exs'] = mass_exs
            maps['mass_gal'] = mass_gal
            
        else:
            print('Note: map_star_rho_insitu and map_star_rho_exsitu not in file. Stellar mass maps skipped.')

            
        # Get the stellar age map if there
        if 'map_star_age_insitu_{}'.format(proj) in self.hdf5_values and 'map_star_age_exsitu_{}'.format(proj) in self.hdf5_values:
            age_ins = self.data['map_star_age_insitu_{}'.format(proj)][idx]
            age_exs = self.data['map_star_age_exsitu_{}'.format(proj)][idx]
            age_gal = (age_ins * mass_ins + age_exs * mass_exs) / (mass_ins + mass_exs)
            age_ins[age_ins == 0.] = np.nan
            age_exs[age_ins == 0.] = np.nan

            if verbose:
                print("\n# (Age_ins/Gyr): {:6.2f}".format(np.nanmedian(age_ins)))
                print("# (Age_exs/Gyr): {:6.2f}".format(np.nanmedian(age_exs)))
            
            maps['age_ins'] = age_ins
            maps['age_exs'] = age_exs
            maps['age_gal'] = age_gal
        else:
            print('map_star_age_insitu and map_star_age_exsitu not in file. Stellar age maps skipped.')
                                                                                   

        # Get the stellar metallicity map if there
        if 'map_star_metallicity_insitu_{}'.format(proj) in self.hdf5_values and 'map_star_metallicity_exsitu_{}'.format(proj) in self.hdf5_values:
            met_ins = self.data['map_star_metallicity_insitu_{}'.format(proj)][idx]
            met_exs = self.data['map_star_metallicity_exsitu_{}'.format(proj)][idx]
            met_gal = (met_ins * mass_ins + met_exs * mass_exs) / (mass_ins + mass_exs)
            met_ins[met_ins == 0.] = np.nan
            met_exs[met_ins == 0.] = np.nan

            if verbose:
                print("# log(Z_ins/Z_sun): {:6.2f}".format(
                    np.log10(np.nanmedian(met_ins / Z_SUN))))
                print("# log(Z_exs/Z_sun): {:6.2f}".format(
                    np.log10(np.nanmedian(met_exs / Z_SUN))))

            maps['met_ins'] = met_ins
            maps['met_exs'] = met_exs
            maps['met_gal'] = met_gal
        else:
            print('map_star_metallicity_insitu and map_star_metallicity_exsitu not in file. Metallicity maps skipped.')
                
        if maps_only:
            return maps

        return info, maps
