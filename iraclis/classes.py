
__all__ = ['variables', 'calibrations','tools', 'DataSet', 'PipelineCounter']

import os
import sys
import glob
import time
import numpy as np
import shutil
import astropy.io.fits as pf
import datetime
import matplotlib.pyplot as plt
import pylightcurve as plc

from scipy.optimize import curve_fit

from iraclis.errors import *
from iraclis.__databases__ import iraclis_data


class Variable:

    def __init__(self, name, keyword=None, value=None, instance=None, kind='p', hdu_image=False, auto_fill=False):

        self.name = name
        if keyword is None:
            keyword = name
        self.keyword = keyword
        self.default_value = value
        self.value = value
        self.instance = instance
        self.kind = kind
        self.hdu_image = hdu_image
        self.auto_fill = auto_fill

    def set(self, value):

        if value is not None:

            if self.instance is None:
                self.value = value

            elif self.hdu_image:
                self.value = value

            elif self.auto_fill:
                if isinstance(value, str):
                    if value in ['auto', 'default', 'default_high', 'default_low', 'default_vlow']:
                        self.value = value
                    else:
                        try:
                            self.value = self.instance(value)
                        except ValueError:
                            print(value)
                            raise IraclisInputError('Input {0} parameter is not valid, '
                                                    '{1} is expected.'.format(self.name, self.instance))
                else:
                    try:
                        self.value = self.instance(value)
                    except ValueError:
                        raise IraclisInputError('Input {0} parameter is not valid, '
                                                '{1} is expected.'.format(self.name, self.instance))
            elif self.instance is bool:

                if isinstance(value, bool):
                    self.value = value
                elif isinstance(value, str):
                    if value == 'True':
                        self.value = True
                    elif value == 'False':
                        self.value = False
                else:
                    raise IraclisInputError('Input {0} parameter is not valid, '
                                            '{1} is expected.'.format(self.name, self.instance))
            else:
                try:
                    self.value = self.instance(value)
                except ValueError:
                        raise IraclisInputError('Input {0} parameter is not valid, '
                                                '{1} is expected.'.format(self.name, self.instance))

    def reset(self):

        self.set(self.default_value)

    def from_fits(self, fits, position=0):

        if self.hdu_image:
            self.set(fits[self.keyword].data)
        else:
            self.set(fits[position].header[self.keyword])

    def to_fits(self, fits, position=0, value=None):

        if value is not None:
            self.set(value)

        if self.hdu_image:
            try:
                fits[self.keyword] = pf.ImageHDU(name=self.keyword, header=fits[self.keyword].header, data=self.value)
            except KeyError:
                fits.append(pf.ImageHDU(name=self.keyword, data=self.value))
        else:
            fits[position].header.set(self.keyword, self.value)

    def from_dictionary(self, dictionary, sub_dictionary=None):

        if isinstance(dictionary, Variable):
            if sub_dictionary is None:
                self.set(dictionary.value[self.keyword])
            else:
                self.set(dictionary.value[sub_dictionary][self.keyword])
        else:
            if sub_dictionary is None:
                self.set(dictionary[self.keyword])
            else:
                self.set(dictionary[sub_dictionary][self.keyword])

    def to_dictionary(self, dictionary, sub_dictionary=None, value=None):

        if value is not None:
            self.set(value)

        if isinstance(dictionary, Variable):
            if sub_dictionary is None:
                dictionary.value[self.keyword] = self.value
            else:
                dictionary.value[sub_dictionary][self.keyword] = self.value
        else:
            if sub_dictionary is None:
                dictionary[self.keyword] = self.value
            else:
                dictionary[sub_dictionary][self.keyword] = self.value

    def custom(self, value=None):

        if value is None:
            return Variable(self.name, self.keyword, self.value, self.instance, self.kind,
                            self.hdu_image, self.auto_fill)
        else:
            return Variable(self.name, self.keyword, value, self.instance, self.kind, self.hdu_image,
                            self.auto_fill)

    def custom_from_fits(self, fits, position=0):

        x = self.custom()
        x.from_fits(fits, position)

        return x

    def custom_from_dictionary(self, dictionary, sub_dictionary=None):

        x = self.custom()
        x.from_dictionary(dictionary, sub_dictionary)

        return x


class Variables:

    def __init__(self):

        # process visit
        self.data_directory = Variable(
            'data_directory', value='.', instance=str, kind='u')
        self.reduction = Variable(
            'reduction', value=True, instance=bool, kind='u')
        self.splitting = Variable(
            'splitting', value=True, instance=bool, kind='u')
        self.extraction = Variable(
            'extraction', value=False, instance=bool, kind='u')
        self.splitting_extraction = Variable(
            'splitting_extraction', value=True, instance=bool, kind='u')
        self.fitting_white = Variable(
            'fitting_white', value=True, instance=bool, kind='u')
        self.fitting_spectrum = Variable(
            'fitting_spectrum', value=True, instance=bool, kind='u')
        self.raw_data_directory = Variable(
            'raw_data_directory', value='raw_data', instance=str)
        self.reduced_data_directory = Variable(
            'reduced_data_directory', value='reduced_data', instance=str)
        self.splitted_data_directory = Variable(
            'splitted_data_directory', value='splitted_data', instance=str)
        self.figures_directory = Variable(
            'figures_directory', value='figures', instance=str)
        self.output_directory = Variable(
            'output_directory', value='results', instance=str)
        self.output_directory_copy = Variable(
            'output_directory_copy', value='False', instance=str, kind='u')
        self.light_curve_file = Variable(
            'light_curve_file', value='extracted_light_curves', instance=str)
        self.fitting_file = Variable(
            'fitting_file', value='fitting_results', instance=str)
        self.apply_bias = Variable(
            'apply_bias', value=True, instance=bool)
        self.apply_linearity = Variable(
            'apply_linearity', value=True, instance=bool)
        self.apply_dark = Variable(
            'apply_dark', value=True, instance=bool)
        self.apply_gain = Variable(
            'apply_gain', value=True, instance=bool)
        self.apply_sky = Variable(
            'apply_sky', value=True, instance=bool)
        self.apply_flat = Variable(
            'apply_flat', value=True, instance=bool)
        self.apply_bpcr = Variable(
            'apply_bpcr', value=True, instance=bool)

        # observation
        self.ra_target = Variable(
            'ra_target', keyword='RA_TARG', instance=float)
        self.dec_target = Variable(
            'dec_target', keyword='DEC_TARG', instance=float)
        self.observation_type = Variable(
            'observation_type', keyword='OBSTYPE', instance=str)
        self.grism = Variable(
            'grism', keyword='FILTER', instance=str)
        self.wfc3_aperture = Variable(
            'wfc3_aperture', keyword='APERTURE', instance=str)
        self.postarg1 = Variable(
            'postarg1', keyword='POSTARG1', instance=float)
        self.scan_direction = Variable(
            'scan_direction', keyword='POSTARG2', instance=float)
        self.exposure_start = Variable(
            'exposure_start', keyword='EXPSTART', instance=float)
        self.exposure_end = Variable(
            'exposure_end', keyword='EXPEND', instance=float)
        self.exposure_time = Variable(
            'exposure_time', keyword='EXPTIME', instance=float)
        self.sub_array_type = Variable(
            'sub_array_type', keyword='SUBTYPE', instance=str)
        self.sub_array_size = Variable(
            'sub_array_size', keyword='SUBARRAY', instance=int)
        self.total_samples = Variable(
            'total_samples', keyword='NSAMP', instance=int)
        self.sampling_sequence = Variable(
            'sampling_sequence', keyword='SAMP_SEQ', instance=str)
        self.sample_number = Variable(
            'sample_number', keyword='SAMPNUM', instance=int)
        self.sample_time = Variable(
            'sample_time', keyword='SAMPTIME', instance=float)

        # timing
        self.bjd_tdb = Variable('bjd_tdb', keyword='BJD_TDB', instance=float)

        # bias and zero-read
        self.zero_read = Variable(
            'zero_read', keyword='ZREAD', instance=np.array, hdu_image=True)
        self.zero_read_error = Variable(
            'zero_read_error', keyword='ZREADE', instance=np.array, hdu_image=True)
        self.zero_read_flux = Variable(
            'zero_read_flux', keyword='ZFLUX', instance=np.array, hdu_image=True)
        self.zero_read_flux_error = Variable(
            'zero_read_flux_error', keyword='ZFLUXE', instance=np.array, hdu_image=True)
        self.reference_pixels_level = Variable(
            'reference_pixels_level', keyword='REFLEV', instance=float)

        # sky
        self.sky_detection_limit = Variable(
            'sky_detection_limit', keyword='SKYLIM', value=2.0, instance=float)
        self.sky_background_level = Variable(
            'sky_background_level', keyword='SKYBGLEV', instance=float)
        self.sky_area = Variable(
            'sky_area', keyword='SKYAREA', instance=np.array, hdu_image=True)

        # bpcr
        self.cr_neighbours = Variable(
            'cr_neighbours', keyword='CRXNB', value=6, instance=int)
        self.cr_detection_limit = Variable(
            'cr_detection_limit', keyword='CRLIM', value=3.0, instance=float)
        self.use_bpcr_fast_mode = Variable(
            'use_bpcr_fast_mode', keyword='FBPCR', value=True, instance=bool)
        self.bpcr_map = Variable(
            'bpcr_map', keyword='BPCRMAP', instance=np.array, hdu_image=True)

        # calibration

        self.reference_pixel_x = Variable('reference_pixel_x', keyword='CRPIX1', instance=int)
        self.reference_pixel_y = Variable('reference_pixel_y', keyword='CRPIX2', instance=int)

        self.comparison_index_forward = Variable(
            'comparison_index_forward', keyword='FCOMP', value=0, instance=int)
        self.comparison_index_reverse = Variable(
            'comparison_index_reverse', keyword='RCOMP', value=1, instance=int)
        self.target_x_offset = Variable(
            'target_x_offset', keyword='TARXOFF', value=0, instance=float, kind='u')
        self.target_y_offset = Variable(
            'target_y_offset', keyword='TARYOFF', value=0, instance=float, kind='u')
        self.use_standard_flat = Variable(
            'use_standard_flat', keyword='STDFLAT', value=True, instance=bool)

        self.spectrum_bottom = Variable(
            'spectrum_bottom', keyword='SPCBTTM', instance=int)
        self.spectrum_top = Variable(
            'spectrum_top', keyword='SPCTOP', instance=int)
        self.spectrum_left = Variable(
            'spectrum_left', keyword='SPCLEFT', instance=int)
        self.spectrum_right = Variable(
            'spectrum_right', keyword='SPCRIGHT', instance=int)
        self.spectrum_scan = Variable(
            'spectrum_scan', keyword='SPCSCAN', instance=bool)
        self.spectrum_direction = Variable(
            'spectrum_direction', keyword='SPECDIR', instance=int)
        self.first_spectrum_bottom = Variable(
            'first_spectrum_bottom', keyword='FSPCBTTM', instance=int)
        self.first_spectrum_top = Variable(
            'first_spectrum_top', keyword='FSPCTOP', instance=int)
        self.first_spectrum_scan = Variable(
            'first_spectrum_scan', keyword='FSPCSCAN', instance=bool)
        self.first_spectrum_direction = Variable(
            'first_spectrum_direction', keyword='FSPECDIR', instance=int)

        self.comparison_x_star = Variable(
            'comparison_x_star', keyword='CMPXSTAR', instance=float)
        self.x_star = Variable(
            'x_star', keyword='XSTAR', instance=float)
        self.x_shift = Variable(
            'x_shift', keyword='XSHIFT', instance=float)
        self.x_shift_error = Variable(
            'x_shift_error', keyword='XSHIFTE', instance=float)

        self.comparison_y_star = Variable(
            'comparison_y_star', keyword='CMPYSTAR', instance=float)
        self.y_star = Variable(
            'y_star', keyword='YSTAR', instance=float)
        self.y_shift = Variable(
            'y_shift', keyword='YSHIFT', instance=float)
        self.y_shift_error = Variable(
            'y_shift_error', keyword='YSHIFTE', instance=float)

        self.scan_length = Variable(
            'scan_length', keyword='LEN', instance=float)
        self.scan_length_error = Variable(
            'scan_length_error', keyword='LEN_ERR', instance=float)

        self.wdpt_constant_coefficient_1 = Variable(
            'wdpt_constant_coefficient_1', keyword='CSCOEFF1', instance=float)
        self.wdpt_constant_coefficient_2 = Variable(
            'wdpt_constant_coefficient_2', keyword='CSCOEFF2', instance=float)
        self.wdpt_constant_coefficient_3 = Variable(
            'wdpt_constant_coefficient_3', keyword='CSCOEFF3', instance=float)
        self.wdpt_slope_coefficient_1 = Variable(
            'wdpt_slope_coefficient_1', keyword='SLCOEFF1', instance=float)
        self.wdpt_slope_coefficient_2 = Variable(
            'wdpt_slope_coefficient_2', keyword='SLCOEFF2', instance=float)
        self.wdpt_slope_coefficient_3 = Variable(
            'wdpt_slope_coefficient_3', keyword='SLCOEFF3', instance=float)

        self.scan_frame = Variable(
            'scan_frame', keyword='SCANMAP', instance=np.array, hdu_image=True)
        self.wavelength_frame = Variable(
            'wavelength_frame', keyword='WMAP', instance=np.array, hdu_image=True)
        self.normalised_wavelength_frame = Variable(
            'normalised_wavelength_frame', keyword='NWMAP', instance=np.array, hdu_image=True)

        # extraction
        self.aperture_lower_extend = Variable(
            'aperture_lower_extend', value=-20.0, instance=float, kind='u')
        self.aperture_upper_extend = Variable(
            'aperture_upper_extend', value=20.0, instance=float, kind='u')
        self.extraction_method = Variable(
            'extraction_method', value='gauss', instance=str, kind='u')
        self.extraction_gauss_sigma = Variable(
            'extraction_gauss_sigma', value=47.0, instance=float, kind='u')
        self.white_lower_wavelength = Variable(
            'white_lower_wavelength', value='default', instance=float, kind='u', auto_fill=True)
        self.white_upper_wavelength = Variable(
            'white_upper_wavelength', value='default', instance=float, kind='u', auto_fill=True)
        self.bins_file = Variable(
            'bins_file', value='default_high', instance=str, kind='u')
        self.bins_number = Variable(
            'bins_number', value=0, instance=int)
        self.bjd_tdb_array = Variable(
            'bjd_tdb_array', value=np.array([]), instance=np.array)
        self.spectrum_direction_array = Variable(
            'spectrum_direction_array', value=np.array([]), instance=np.array)
        self.sky_background_level_array = Variable(
            'sky_background_level_array', value=np.array([]), instance=np.array)
        self.x_star_array = Variable(
            'x_star_array', value=np.array([]), instance=np.array)
        self.x_shift_error_array = Variable(
            'x_shift_error_array', value=np.array([]), instance=np.array)
        self.y_star_array = Variable(
            'y_star_array', value=np.array([]), instance=np.array)
        self.y_shift_error_array = Variable(
            'y_shift_error_array', value=np.array([]), instance=np.array)
        self.scan_length_array = Variable(
            'scan_length_array', value=np.array([]), instance=np.array)
        self.scan_length_error_array = Variable(
            'scan_length_error_array', value=np.array([]), instance=np.array)
        self.lower_wavelength = Variable(
            'lower_wavelength', instance=float)
        self.upper_wavelength = Variable(
            'upper_wavelength', instance=float)
        self.ldc1 = Variable(
            'ldc1', instance=float, auto_fill=True)
        self.ldc2 = Variable(
            'ldc2', instance=float, auto_fill=True)
        self.ldc3 = Variable(
            'ldc3', instance=float, auto_fill=True)
        self.ldc4 = Variable(
            'ldc4', instance=float, auto_fill=True)
        self.flux_array = Variable(
            'flux', value=np.array([]), instance=np.array)
        self.error_array = Variable(
            'error', value=np.array([]), instance=np.array)
        self.ph_error_array = Variable(
            'ph_error', value=np.array([]), instance=np.array)
        self.white_dictionary = Variable(
            'white', value={}, instance=dict)
        self.light_curve_split = Variable(
            'light_curve_split', 'split_', instance=str)

        # fitting
        self.apply_up_down_stream_correction = Variable('apply_up_down_stream_correction', 'USDS', True, bool, 'u')
        self.exclude_initial_orbits = Variable('exclude_initial_orbits', 'EXCL1', 1, int, 'u')
        self.exclude_final_orbits = Variable('exclude_final_orbits', 'EXCL2', 0, int, 'u')
        self.exclude_initial_orbit_points = Variable('exclude_initial_orbit_points', 'EXCP1', 0, int, 'u')
        self.white_ldc1 = Variable('white_ldc1', 'WHTLDC1', 'default', float, 'u', auto_fill=True)
        self.white_ldc2 = Variable('white_ldc2', 'WHTLDC2', 'default', float, 'u', auto_fill=True)
        self.white_ldc3 = Variable('white_ldc3', 'WHTLDC3', 'default', float, 'u', auto_fill=True)
        self.white_ldc4 = Variable('white_ldc4', 'WHTLDC4', 'default', float, 'u', auto_fill=True)
        self.fit_ldc1 = Variable('fit_ldc1', 'FITLDC1', False, bool, 'u')
        self.fit_ldc2 = Variable('fit_ldc2', 'FITLDC2', False, bool, 'u')
        self.fit_ldc3 = Variable('fit_ldc3', 'FITLDC3', False, bool, 'u')
        self.fit_ldc4 = Variable('fit_ldc4', 'FITLDC4', False, bool, 'u')
        self.planet = Variable('planet', 'PLANET', 'auto', str, 'u', auto_fill=True)
        self.method = Variable('method', 'METHOD', 'claret', str, 'u')
        self.star_teff = Variable('star_teff', 'TEFF', 'auto', float, 'u', auto_fill=True)
        self.star_logg = Variable('star_logg', 'LOGG', 'auto', float, 'u', auto_fill=True)
        self.star_meta = Variable('star_meta', 'META', 'auto', float, 'u', auto_fill=True)
        self.rp_over_rs = Variable('rp_over_rs ', 'RP', 'auto', float, 'u', auto_fill=True)
        self.fp_over_fs = Variable('fp_over_fs', 'FP', 'auto', float, 'u', auto_fill=True)
        self.period = Variable('period', 'P', 'auto', float, 'u', auto_fill=True)
        self.sma_over_rs = Variable('sma_over_rs', 'A', 'auto', float, 'u', auto_fill=True)
        self.fit_sma_over_rs = Variable('fit_sma_over_rs', 'FITA', False, bool, 'u')
        self.eccentricity = Variable('eccentricity', 'E', 'auto', float, 'u', auto_fill=True)
        self.inclination = Variable('inclination', 'I', 'auto', float, 'u', auto_fill=True)
        self.fit_inclination = Variable('fit_inclination', 'FITI', False, bool, 'u')
        self.periastron = Variable('periastron', 'W', 'auto', float, 'u', auto_fill=True)
        self.mid_time = Variable('mid_time', 'MT', 'auto', float, 'u', auto_fill=True)
        self.fit_mid_time = Variable('fit_mid_time', 'FITMT', True, bool, 'u')
        self.second_order_ramp = Variable('second_order_ramp', 'RAMP2', False, bool, 'u')
        self.first_orbit_ramp = Variable('first_orbit_ramp', 'FOR', True, bool, 'u')
        self.mid_orbit_ramps = Variable('mid_orbit_ramps', 'MOR', True, bool, 'u')
        self.mcmc_iterations = Variable('mcmc_iterations', 'MCMCI', 300000, int, 'u')
        self.mcmc_walkers = Variable('mcmc_walkers', 'MCMCW', 50, int, 'u')
        self.mcmc_burned_iterations = Variable('mcmc_burned_iterations', 'MCMCB', 200000, int, 'u')
        self.spectral_mcmc_iterations = Variable('spectral_mcmc_iterations', 'SPMCMCI', 60000, int, 'u')
        self.spectral_mcmc_walkers = Variable('spectral_mcmc_walkers', 'SPMCMCW', 50, int, 'u')
        self.spectral_mcmc_burned_iterations = Variable('spectral_mcmc_burned_iterations', 'SPMCMCB', 10000, int, 'u')

    def from_parameters_file(self, parameters_file=None, data_directory=None):

        # check user inputs

        variables.reset()

        if parameters_file and data_directory:

            raise IraclisInputError('You should give as input either a parameters file '
                                    '(which includes a data directory path), '
                                    'OR a data directory path '
                                    '(if you want to run with default parameters).')

        elif parameters_file:

            parameters_file = os.path.abspath(parameters_file)
            if not os.path.isfile(parameters_file):
                raise IraclisFileError('No such file: ' + parameters_file)
            print('Loading parameters file: {} ...'.format(os.path.split(parameters_file)[1]))

            parameters_in_file = []

            for i in open(parameters_file).readlines():
                if len(i.split()) > 0:
                    if i.split()[0][0] != '#':
                        if i.split()[0] in vars(self):
                            parameters_in_file.append(i.split()[0])
                            vars(self)[i.split()[0]].set(i.split()[1])

            for i in vars(self):
                if not isinstance(vars(self)[i], list):
                    if i not in parameters_in_file:
                        if vars(self)[i].kind == 'u':
                            print('WARNING: Parameter not in file, setting {0}={1}'.format(
                                i, vars(self)[i].default_value))

        elif data_directory:

            self.data_directory.set(data_directory)

        else:

            raise IraclisInputError('You should give as input either a parameters file '
                                    '(which includes a data directory path), '
                                    'OR a data directory path '
                                    '(if you want to run in default parameters).')

    def overwrite(self, procedure=None, par_string=None):

        if procedure:
            self.reduction.set(bool(int(procedure[0])))
            self.splitting.set(bool(int(procedure[1])))
            self.extraction.set(bool(int(procedure[2])))
            self.splitting_extraction.set(bool(int(procedure[3])))
            self.fitting_white.set(bool(int(procedure[4])))
            self.fitting_spectrum.set(bool(int(procedure[5])))

        if par_string:
            for i in par_string.split(','):
                if len(i.split('=')) > 0:
                    if i.split('=')[0] in vars(self):
                        vars(self)[i.split('=')[0]].set(i.split('=')[1])

    def set_binning(self, input_data, white_lower_wavelength, white_upper_wavelength, white_ldc1,
                    white_ldc2, white_ldc3, white_ldc4, bins_file):

        grism = self.grism.custom()
        sub_array_size = self.sub_array_size.custom()

        if isinstance(input_data, dict):
            grism.from_dictionary(input_data)
            sub_array_size.from_dictionary(input_data)
        elif input_data.splitted:
            grism.from_fits(input_data.spectroscopic_images[0][0])
            sub_array_size.from_fits(input_data.spectroscopic_images[0][0])
        else:
            grism.from_fits(input_data.spectroscopic_images[0])
            sub_array_size.from_fits(input_data.spectroscopic_images[0])

        if grism.value == 'G141':
            if sub_array_size.value < 200:
                if white_lower_wavelength == 'default':
                    white_lower_wavelength = 11000
                if white_upper_wavelength == 'default':
                    white_upper_wavelength = 16500
            else:
                if white_lower_wavelength == 'default':
                    white_lower_wavelength = 10800
                if white_upper_wavelength == 'default':
                    white_upper_wavelength = 16800

        elif grism.value == 'G102':
            if white_lower_wavelength == 'default':
                white_lower_wavelength = 8000
            if white_upper_wavelength == 'default':
                white_upper_wavelength = 11500

        if os.path.isfile(os.path.abspath(bins_file)):
            bins_file = os.path.abspath(bins_file)
            print('Loading bins file: {} ...'.format(os.path.split(bins_file)[1]))

            if len(np.loadtxt(bins_file, unpack=True)) == 6:
                bins_lower_wavelength, bins_upper_wavelength, bins_ldc1, bins_ldc2, bins_ldc3, bins_ldc4 = \
                    np.loadtxt(bins_file, unpack=True)
                bins_ldc1 = list(bins_ldc1)
                bins_ldc2 = list(bins_ldc2)
                bins_ldc3 = list(bins_ldc3)
                bins_ldc4 = list(bins_ldc4)

            else:
                raise IraclisFileError('Not valid bins file.')

        elif bins_file == 'default_high':
            print('Loading bins file: {} ...'.format('default_high.txt'))
            if grism.value == 'G141':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_high.txt')
            elif grism.value == 'G102':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_high_g102.txt')

            bins_lower_wavelength, bins_upper_wavelength = np.loadtxt(bins_file, unpack=True)
            bins_ldc1 = ['default_high' for ff in bins_lower_wavelength]
            bins_ldc2 = ['default_high' for ff in bins_lower_wavelength]
            bins_ldc3 = ['default_high' for ff in bins_lower_wavelength]
            bins_ldc4 = ['default_high' for ff in bins_lower_wavelength]

        elif bins_file == 'default_low':
            print('Loading bins file: {} ...'.format('default_low.txt'))
            if grism.value == 'G141':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_low.txt')
            elif grism.value == 'G102':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_low_g102.txt')

            bins_lower_wavelength, bins_upper_wavelength = np.loadtxt(bins_file, unpack=True)
            bins_ldc1 = ['default_low' for ff in bins_lower_wavelength]
            bins_ldc2 = ['default_low' for ff in bins_lower_wavelength]
            bins_ldc3 = ['default_low' for ff in bins_lower_wavelength]
            bins_ldc4 = ['default_low' for ff in bins_lower_wavelength]

        elif bins_file == 'default_vlow':
            print('Loading bins file: {} ...'.format('default_vlow.txt'))
            if grism.value == 'G141':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_vlow.txt')
            elif grism.value == 'G102':
                bins_file = os.path.join(iraclis_data.hstwfc3(), 'default_vlow_g102.txt')

            bins_lower_wavelength, bins_upper_wavelength = np.loadtxt(bins_file, unpack=True)
            bins_ldc1 = ['default_vlow' for ff in bins_lower_wavelength]
            bins_ldc2 = ['default_vlow' for ff in bins_lower_wavelength]
            bins_ldc3 = ['default_vlow' for ff in bins_lower_wavelength]
            bins_ldc4 = ['default_vlow' for ff in bins_lower_wavelength]

        else:
            raise IraclisFileError('No bins file found')

        lower_wavelength = self.lower_wavelength.custom()
        upper_wavelength = self.upper_wavelength.custom()
        flux_array = self.flux_array.custom()
        error_array = self.error_array.custom()
        ph_error_array = self.ph_error_array.custom()
        ldc1 = self.ldc1.custom()
        ldc2 = self.ldc2.custom()
        ldc3 = self.ldc3.custom()
        ldc4 = self.ldc4.custom()

        white_dictionary = Variable('white', value={}, instance=dict)
        bins_dictionaries = [Variable('bin_{0}'.format(str(ff + 1).zfill(2)), value={}, instance=dict)
                             for ff in range(len(bins_lower_wavelength))]

        for i, j in enumerate([white_dictionary] + bins_dictionaries):
            flux_array.to_dictionary(j)
            error_array.to_dictionary(j)
            ph_error_array.to_dictionary(j)
            lower_wavelength.to_dictionary(j, value=([white_lower_wavelength] + list(bins_lower_wavelength))[i])
            upper_wavelength.to_dictionary(j, value=([white_upper_wavelength] + list(bins_upper_wavelength))[i])
            ldc1.to_dictionary(j, value=([white_ldc1] + bins_ldc1)[i])
            ldc2.to_dictionary(j, value=([white_ldc2] + bins_ldc2)[i])
            ldc3.to_dictionary(j, value=([white_ldc3] + bins_ldc3)[i])
            ldc4.to_dictionary(j, value=([white_ldc4] + bins_ldc4)[i])

        return white_dictionary, bins_dictionaries

    def reset(self):

        for i in vars(self):
            if isinstance(vars(self)[i], list):
                for j in vars(self)[i]:
                    j.reset()
            else:
                vars(self)[i].reset()

    def save(self, export_file):

        w = open(export_file, 'w')

        for i in vars(self):
            if not isinstance(vars(self)[i], list):
                if vars(self)[i].kind == 'u':
                    w.write('{0}{1}{2}\n'.format(i, ' ' * (35 - len(i)), vars(self)[i].value))

        w.close()


variables = Variables()


# calibration variables


class Calibration:

    def __init__(self, calibration_files_id, keyword, value=None, matching_keywords=None,
                 hdu_image=False, dark_mode=False):

        self.calibration_files_id = calibration_files_id
        self.calibration_data = None
        self.keyword = keyword
        self.value = value
        self.hdu_image = hdu_image
        self.dark_mode = dark_mode
        self.matching_keywords = matching_keywords

    def match(self, input_data):

        if isinstance(input_data, pf.HDUList):
            pass
        elif input_data.splitted:
            input_data = input_data.spectroscopic_images[0][0]
        else:
            input_data = input_data.spectroscopic_images[0]

        calibration_directory_path = iraclis_data.hstwfc3()

        if not self.calibration_data:

            self.calibration_data = \
                [pf.open(pp) for pp in
                 glob.glob(os.path.join(calibration_directory_path, '*{0}*.fits'.format(self.calibration_files_id)))]

        if len(self.calibration_data) == 1 or not self.matching_keywords:

            test = True
            selected_calibration_data = self.calibration_data[0]

        else:

            test = False
            selected_calibration_data = self.calibration_data[0]
            for selected_calibration_data in self.calibration_data:
                if ([input_data[0].header[pp] for pp in self.matching_keywords] ==
                        [selected_calibration_data[0].header[pp] for pp in self.matching_keywords]):
                    test = True
                    break

        if test:

            if self.hdu_image:

                if self.dark_mode:

                    dark_list = np.array([selected_calibration_data[0]])
                    sample_numbers = [input_data[pp].header[variables.sample_number.keyword]
                                      for pp in plc.fits_sci(input_data)]

                    for i in plc.fits_sci(selected_calibration_data):
                        if (selected_calibration_data[i].header[
                                variables.sample_number.keyword] in sample_numbers
                            or selected_calibration_data[i].header[
                                variables.sample_number.keyword] == min(sample_numbers) - 1):

                            dark_list = np.append(dark_list, selected_calibration_data[i])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 1])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 2])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 3])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 4])

                    return pf.HDUList(list(dark_list))

                else:

                    array = (np.ones_like(selected_calibration_data[self.keyword].data) *
                             selected_calibration_data[self.keyword].data)

                    crop1 = int(len(np.array(array)) / 2 - len(input_data[plc.fits_sci(input_data)[0]].data) / 2)
                    crop2 = int(len(np.array(array)) / 2 + len(input_data[plc.fits_sci(input_data)[0]].data) / 2)

                    return array[crop1:crop2, crop1:crop2]

            else:
                return (np.ones_like(selected_calibration_data[0].header[self.keyword]) *
                        selected_calibration_data[0].header[self.keyword])

        else:
            raise IraclisLibraryError('No suitable dark current calibration file found, \n'
                                      'you must update your library manually. Search at: \n\n'
                                      '\thttps://hst-crds.stsci.edu/ \n\n'
                                      'for the latest dark current calibration file for\n\n'
                                      '\tSAMP_SEQ: {0} \n\tSUBTYTE: {1}'.format(input_data[0].header['SAMP_SEQ'],
                                                                                input_data[0].header['SUBTYPE']))


class Calibrations:

    def __init__(self):

        self.super_zero_read = Calibration(
            'lin_custom', 'ZSCI', hdu_image=True)
        self.super_zero_read_error = Calibration(
            'lin_custom', 'ZERR', hdu_image=True)
        self.ccd_gain = Calibration(
            'ccd_custom', 'CCDGAIN', hdu_image=True)
        self.ccd_read_noise = Calibration(
            'ccd_custom', 'CCDREADN', hdu_image=True)
        self.linearity_coefficient_1 = Calibration(
            'lin_custom', 'LINC1', hdu_image=True)
        self.linearity_coefficient_2 = Calibration(
            'lin_custom', 'LINC2', hdu_image=True)
        self.linearity_coefficient_3 = Calibration(
            'lin_custom', 'LINC3', hdu_image=True)
        self.linearity_coefficient_4 = Calibration(
            'lin_custom', 'LINC4', hdu_image=True)
        self.super_dark = Calibration(
            'drk', 'SCI', dark_mode=True, hdu_image=True,
            matching_keywords=[variables.sampling_sequence.keyword, variables.sub_array_type.keyword])
        self.mean_gain = Calibration(
            'ccd_custom', 'MGAIN')
        self.gain_profile = Calibration(
            'pfl', 'SCI', hdu_image=True)
        self.master_sky = Calibration(
            'sky', 'PRIMARY', hdu_image=True, matching_keywords=[variables.grism.keyword])
        self.flat_field_coefficient_1 = Calibration(
            'flat', 'FFC1', hdu_image=True, matching_keywords=[variables.grism.keyword])
        self.flat_field_coefficient_2 = Calibration(
            'flat', 'FFC2', hdu_image=True, matching_keywords=[variables.grism.keyword])
        self.flat_field_coefficient_3 = Calibration(
            'flat', 'FFC3', hdu_image=True, matching_keywords=[variables.grism.keyword])
        self.flat_field_coefficient_4 = Calibration(
            'flat', 'FFC4', hdu_image=True, matching_keywords=[variables.grism.keyword])
        self.flat_field_min_wavelength = Calibration(
            'flat', 'WMIN', matching_keywords=[variables.grism.keyword])
        self.flat_field_max_wavelength = Calibration(
            'flat', 'WMAX', matching_keywords=[variables.grism.keyword])
        self.bad_pixels = Calibration(
            'bpx_custom', 'PRIMARY', hdu_image=True)
        self.trace_at0 = Calibration(
            'calib', 'AT0', matching_keywords=[variables.grism.keyword])
        self.trace_at1 = Calibration(
            'calib', 'AT1', matching_keywords=[variables.grism.keyword])
        self.trace_at2 = Calibration(
            'calib', 'AT2', matching_keywords=[variables.grism.keyword])
        self.trace_at3 = Calibration(
            'calib', 'AT3', matching_keywords=[variables.grism.keyword])
        self.trace_at4 = Calibration(
            'calib', 'AT4', matching_keywords=[variables.grism.keyword])
        self.trace_at5 = Calibration(
            'calib', 'AT5', matching_keywords=[variables.grism.keyword])
        self.trace_bt0 = Calibration(
            'calib', 'BT0', matching_keywords=[variables.grism.keyword])
        self.trace_bt1 = Calibration(
            'calib', 'BT1', matching_keywords=[variables.grism.keyword])
        self.trace_bt2 = Calibration(
            'calib', 'BT2', matching_keywords=[variables.grism.keyword])
        self.wsol_aw0 = Calibration(
            'calib', 'AW0', matching_keywords=[variables.grism.keyword])
        self.wsol_aw1 = Calibration(
            'calib', 'AW1', matching_keywords=[variables.grism.keyword])
        self.wsol_aw2 = Calibration(
            'calib', 'AW2', matching_keywords=[variables.grism.keyword])
        self.wsol_aw3 = Calibration(
            'calib', 'AW3', matching_keywords=[variables.grism.keyword])
        self.wsol_aw4 = Calibration(
            'calib', 'AW4', matching_keywords=[variables.grism.keyword])
        self.wsol_aw5 = Calibration(
            'calib', 'AW5', matching_keywords=[variables.grism.keyword])
        self.wsol_bw0 = Calibration(
            'calib', 'BW0', matching_keywords=[variables.grism.keyword])
        self.wsol_bw1 = Calibration(
            'calib', 'BW1', matching_keywords=[variables.grism.keyword])
        self.wsol_bw2 = Calibration(
            'calib', 'BW2', matching_keywords=[variables.grism.keyword])


calibrations = Calibrations()


class DataSet:

    def __init__(self, input_data, direct_image=None):

        if direct_image is not None:
            if isinstance(direct_image, str):
                if os.path.isfile(direct_image):
                    direct_image = pf.open(direct_image)
                else:
                    raise IraclisFileError('No such file {0}'.format(input_data))
            else:
                raise IraclisFileError('Please give a file name for the direct image or leave it empty.')

        if not isinstance(input_data, str):
            raise IraclisFileError('Please give a file or directory name for the input data.')

        elif input_data == 'empty':
            self.file_names = []
            self.spectroscopic_images = []
            self.direct_image = [1]
            self.splitted = False
            self._data_set_directory_path = None

        elif os.path.isfile(input_data):

            self.file_names = []
            self.spectroscopic_images = [pf.open(input_data)]
            if not direct_image:
                self.direct_image = None
            else:
                self.direct_image = direct_image
            self.splitted = False
            self._data_set_directory_path = None

            if self.spectroscopic_images[0][0].header[variables.observation_type.keyword] != 'SPECTROSCOPIC':
                raise IraclisFileError('A single direct image is not s valid input dataset.')

        elif os.path.isdir(input_data) and len(glob.glob(os.path.join(input_data, '*', ''))) == 0:

            nsamp = []
            final_list = []
            direct_image = None

            files = sorted(glob.glob(os.path.join(input_data, '*.fits')))
            for i in files:
                with pf.open(i) as j:
                    if j[0].header[variables.observation_type.keyword] == 'SPECTROSCOPIC':
                        final_list.append([j[0].header[variables.exposure_start.keyword],
                                           os.path.split(i)[1], plc.copy_fits(j)])
                        nsamp.append(j[0].header[variables.total_samples.keyword])
                    elif not direct_image:
                        direct_image = plc.copy_fits(j[:2])

            nsamps = [int(np.median(np.array(nsamp)))]

            final_list.sort()
            list_of_times, list_of_files, list_of_fits = np.swapaxes(np.array(final_list, dtype=object), 0, 1)

            outliers = True
            while outliers:
                outliers = False
                for i in range(len(list_of_fits)):
                    if list_of_fits[i][0].header[variables.total_samples.keyword] not in nsamps:
                        list_of_fits = np.delete(list_of_fits, i)
                        list_of_files = np.delete(list_of_files, i)
                        outliers = True
                        break

            self.file_names = list_of_files
            self.spectroscopic_images = list_of_fits
            self.direct_image = direct_image
            self.splitted = False
            self._data_set_directory_path = input_data

        elif os.path.isdir(input_data) and len(glob.glob(os.path.join(input_data, '*', ''))) > 0:

            self.file_names = []
            self.spectroscopic_images = []
            self.direct_image = []
            self.splitted = True
            self._data_set_directory_path = input_data

            for input_data in sorted(glob.glob(os.path.join(input_data, '*', ''))):

                nsamp = []
                final_list = []
                direct_image = None

                files = sorted(glob.glob(os.path.join(input_data, '*.fits')))
                for i in files:
                    with pf.open(i) as j:
                        if j[0].header[variables.observation_type.keyword] == 'SPECTROSCOPIC':
                            final_list.append([j[0].header[variables.exposure_start.keyword],
                                               os.path.split(i)[1], plc.copy_fits(j)])
                            nsamp.append(j[0].header[variables.total_samples.keyword])
                        elif not direct_image:
                            direct_image = plc.copy_fits(j[:2])

                nsamps = [int(np.median(np.array(nsamp)))]

                final_list.sort()
                list_of_times, list_of_files, list_of_fits = np.swapaxes(final_list, 0, 1)

                outliers = True
                while outliers:
                    outliers = False
                    for i in range(len(list_of_fits)):
                        if list_of_fits[i][0].header[variables.total_samples.keyword] not in nsamps:
                            list_of_fits = np.delete(list_of_fits, i)
                            list_of_files = np.delete(list_of_files, i)
                            outliers = True
                            break

                self.file_names = list_of_files
                self.spectroscopic_images.append(list_of_fits)
                self.direct_image = direct_image

        else:
            raise IraclisFileError('No such file or directory: {0}'.format(input_data))

        if not self.direct_image:
            raise IraclisFileError('A direct image is necessary.')

    def save(self, export_directory, arrange=True):

        if os.path.isdir(export_directory):
            backup = '{0}_{1}'.format(export_directory, time.strftime('%y-%m-%d_%H-%M-%S'))
            shutil.copytree(export_directory, backup)
            shutil.rmtree(export_directory)

        os.mkdir(export_directory)

        if arrange:
            for i in range(len(self.file_names)):
                date = str(self.spectroscopic_images[i][0].header['DATE-OBS'])
                obs_time = str(self.spectroscopic_images[i][0].header['TIME-OBS'])
                obs_time = '-'.join(obs_time.split(':'))
                if self.file_names[i].split('_')[0] != date or self.file_names[i].split('_')[1] != obs_time:
                    self.file_names[i] = '{0}_{1}_{2}'.format(date, obs_time, os.path.split(self.file_names[i])[1])

        for i in range(len(self.file_names)):

            copy_of_file = plc.copy_fits(self.spectroscopic_images[i])
            try:
                copy_of_file.writeto(os.path.join(export_directory, self.file_names[i]), output_verify='fix')
            except pf.VerifyError as e:
                hdu = int(str(e.args)[4:-4].split('\\n')[1].replace('HDU ', '').replace(':', ''))
                card = int(str(e.args)[4:-4].split('\\n')[2].replace('Card ', '').replace(':', ''))
                del copy_of_file[hdu].header[card]
                copy_of_file.writeto(os.path.join(export_directory, self.file_names[i]), output_verify='fix')

        copy_of_file = plc.copy_fits(self.direct_image)
        try:
            copy_of_file.writeto(os.path.join(export_directory, 'direct_image.fits'), output_verify='fix')
        except pf.VerifyError as e:
            hdu = int(str(e.args)[4:-4].split('\\n')[1].replace('HDU ', '').replace(':', ''))
            card = int(str(e.args)[4:-4].split('\\n')[2].replace('Card ', '').replace(':', ''))
            del copy_of_file[hdu].header[card]
            copy_of_file.writeto(os.path.join(export_directory, 'direct_image.fits'), output_verify='fix')

    def copy_split(self, split_number):

        x = DataSet('empty')
        x.spectroscopic_images = self.spectroscopic_images[split_number]

        return x


# pipeline counter


class PipelineCounter:

    def __init__(self, task, total_iterations, show_every=1):

        self.task = '{0}{1}'.format(task, '.' * (15 - len(task)))
        self.current_iteration = 0
        self.total_iterations = int(total_iterations)
        self.start_time = time.time()
        self.show = 0
        self.show_every = int(show_every)

        if self.total_iterations == 1:
            self.show_every = 10

    def update(self):

        self.current_iteration += 1
        self.show += 1.0 / self.show_every

        out_of = '{0}{1} / {2}'.format(' ' * (len(str(self.total_iterations)) - len(str(self.current_iteration))),
                                       str(self.current_iteration), str(self.total_iterations))

        delta_time = time.time() - self.start_time

        time_left = str(datetime.timedelta(
            seconds=int((self.total_iterations - self.current_iteration) * delta_time / self.current_iteration)))

        total_time = str(datetime.timedelta(seconds=int(delta_time)))

        if int(self.show):
            sys.stdout.write('\r\033[K')
            sys.stdout.write('{0}: {1}   time left: {2}   total time: {3}'.format(
                self.task, out_of, time_left, total_time))
            sys.stdout.flush()
            self.show = 0

        if self.current_iteration == self.total_iterations and self.total_iterations > 1:
            print('')


class Tools:

    def __init__(self):
        pass

    def central_crop(self, original_array, destination_fits):

        crop1 = len(original_array) / 2 - len(destination_fits[1].data) / 2
        crop2 = len(original_array) / 2 + len(destination_fits[1].data) / 2

        return original_array[crop1:crop2, crop1:crop2]

    def box(self, x_arr, x0, aa, bb, cc):
        return aa * np.exp(-np.power(np.power(x0 - x_arr, 2) / (2 * (bb ** 2)), cc))

    def fit_box(self, datax, datay):

        minlim = datax[np.argmax(datay[5:] - datay[:-5])]
        maxlim = datax[np.argmin(datay[5:] - datay[:-5])]

        for expand in range(10, 30):
            datax = np.append(datax, max(datax) + expand)
            datay = np.append(datay, 0)
            datax = np.append(datax, min(datax) - expand)
            datay = np.append(datay, 0)

        popt, pcov = curve_fit(self.box, datax, datay,
                               p0=[0.5 * (maxlim + minlim), max(datay), 0.5 * (maxlim - minlim), 1.0])

        center = popt[0]

        center_err = np.sqrt(pcov[0][0])

        fwhm = popt[2] * 2.0 * np.sqrt(2.0 * (np.log(2) ** (1.0 / popt[3])))

        s, c = popt[2], popt[3]
        ss, cs = np.sqrt(pcov[2][2]), np.sqrt(pcov[3][3])
        fwhm_err = np.sqrt(8.0 * (ss ** 2) * (np.log(2.0) ** (1.0 / c))
                           + (2.0 * (cs ** 2) * (s ** 2) * (np.log(2.0) ** (1.0 / c))
                              * (np.log(np.log(2.0)) ** 2)) / (c ** 4))

        return center, center_err, fwhm, fwhm_err, popt

    def fit_line(self, datax, datay):

        mx = np.mean(datax)
        my = np.mean(datay)

        ssxx = np.sum((datax - mx) ** 2)
        ssyy = np.sum((datay - my) ** 2)
        ssxy = np.sum((datax - mx) * (datay - my))

        bb = ssxy / ssxx
        aa = my - bb * mx

        n = len(datax)
        sss = np.sqrt((ssyy - bb * ssxy) / (n - 2))

        aerr = sss * np.sqrt(1.0 / n + (mx ** 2) / ssxx)
        berr = sss / np.sqrt(ssxx)

        return aa, bb, aerr, berr

    def correlation(self, x, y):
        n = len(x)
        mx = np.mean(x)
        sx = np.std(x)
        my = np.mean(y)
        sy = np.std(y)
        return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)

    def save_figure(self, directory, name, figure=None):

        if not os.path.isdir(directory):
            os.makedirs(directory)

        if not figure:
            plt.savefig(directory + '/{0}.pdf'.format(name), bbox_inches='tight')
            plt.close()
        else:
            figure.savefig(directory + '/{0}.pdf'.format(name), bbox_inches='tight')


    def adjust_ticks(self):
        xlim1, xlim2 = plt.xlim()
        xticks = plt.xticks()[0]
        digit = max([len(ff) for ff in str(np.abs(xticks) - np.int_(np.abs(xticks)))[1:-1].split()]) - 2
        xticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in xticks]
        dxticks = xticks[1] - xticks[0]

        if xticks[1] - dxticks / 2 < xlim1:
            new_xlim1 = xticks[1] - dxticks / 2
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]
        else:
            new_xlim1 = xticks[0] - dxticks / 2

        if xticks[-2] + dxticks / 2 > xlim2:
            new_xlim2 = xticks[-2] + dxticks / 2
            xticks = xticks[:-1]
            xticklabels = xticklabels[:-1]
        else:
            new_xlim2 = xticks[-1] + dxticks / 2

        plt.xticks(xticks, xticklabels, fontsize=15)
        plt.xlim(new_xlim1, new_xlim2)

        ylim1, ylim2 = plt.ylim()
        yticks = plt.yticks()[0]
        digit = max([len(ff) for ff in str(np.abs(yticks) - np.int_(np.abs(yticks)))[1:-1].split()]) - 2
        yticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in yticks]
        dyticks = yticks[1] - yticks[0]

        if yticks[1] - dyticks / 2 < ylim1:
            new_ylim1 = yticks[1] - dyticks / 2
            yticks = yticks[1:]
            yticklabels = yticklabels[1:]
        else:
            new_ylim1 = yticks[0] - dyticks / 2

        if yticks[-2] + dyticks / 2 > ylim2:
            new_ylim2 = yticks[-2] + dyticks / 2
            yticks = yticks[:-1]
            yticklabels = yticklabels[:-1]
        else:
            new_ylim2 = yticks[-1] + dyticks / 2

        plt.yticks(yticks, yticklabels, fontsize=15)
        plt.ylim(new_ylim1, new_ylim2)


    def adjust_ticks_ax(self, ax):
        xlim1, xlim2 = ax.get_xlim()
        xticks = ax.get_xticks()
        digit = max([len(ff) for ff in str(np.abs(xticks) - np.int_(np.abs(xticks)))[1:-1].split()]) - 2
        xticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in xticks]
        dxticks = xticks[1] - xticks[0]

        if xticks[1] - dxticks / 2 < xlim1:
            new_xlim1 = xticks[1] - dxticks / 2
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]
        else:
            new_xlim1 = xticks[0] - dxticks / 2

        if xticks[-2] + dxticks / 2 > xlim2:
            new_xlim2 = xticks[-2] + dxticks / 2
            xticks = xticks[:-1]
            xticklabels = xticklabels[:-1]
        else:
            new_xlim2 = xticks[-1] + dxticks / 2

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=15)
        ax.set_xlim(new_xlim1, new_xlim2)

        ylim1, ylim2 = ax.get_ylim()
        yticks = ax.get_yticks()
        digit = max([len(ff) for ff in str(np.abs(yticks) - np.int_(np.abs(yticks)))[1:-1].split()]) - 2
        yticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in yticks]
        dyticks = yticks[1] - yticks[0]

        if yticks[1] - dyticks / 2 < ylim1:
            new_ylim1 = yticks[1] - dyticks / 2
            yticks = yticks[1:]
            yticklabels = yticklabels[1:]
        else:
            new_ylim1 = yticks[0] - dyticks / 2

        if yticks[-2] + dyticks / 2 > ylim2:
            new_ylim2 = yticks[-2] + dyticks / 2
            yticks = yticks[:-1]
            yticklabels = yticklabels[:-1]
        else:
            new_ylim2 = yticks[-1] + dyticks / 2

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=15)
        ax.set_ylim(new_ylim1, new_ylim2)

tools = Tools()