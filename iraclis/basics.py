import warnings
warnings.filterwarnings("ignore",
                        message='Matplotlib is building the font cache using fc-list. This may take a moment.')
warnings.filterwarnings("ignore",
                        message='The installed version of numexpr 2.4.4 is not supported in pandas and will be not be used')
warnings.filterwarnings("ignore",
                        message='PyFITS is deprecated, please use astropy.io.fits')

import glob
import os
import urllib
import shutil
import time
import sys
import datetime
import cPickle as pickle
import gzip
import socket

import numpy as np

import scipy
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate

import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
try:
    plt.plot([0], [0], 'o')
    plt.close()
except:
    plt.switch_backend('agg')
import matplotlib.patches as patches

from astropy.io import fits as pf
import ephem
import pylightcurve as plc

import tools


class PYWFC3Error(BaseException):
    pass


class PYWFC3LibraryError(PYWFC3Error):
    pass


class PYWFC3FileError(PYWFC3Error):
    pass


class PYWFC3ProcessError(PYWFC3Error):
    pass


class PYWFC3InputError(PYWFC3Error):
    pass


# useful functions

def open_lightcurve(lc_file):
    return pickle.load(open(lc_file))


def save_lightcurve(lightcurve, export_lc_file):
    pickle.dump(lightcurve, open(export_lc_file, "wb"))


def sci(fits):

    return np.where(np.array([ff.name for ff in fits]) == 'SCI')[0]


def err(fits):

    return np.where(np.array([ff.name for ff in fits]) == 'ERR')[0]


def sci_err(fits):

    return np.swapaxes([sci(fits), err(fits)], 0, 1)


def fits_list(input_data):
    """
    Input data validation.

    Validates the input_data object which must be a HDUList or a DataSet.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    validated_fits_list : list of HDUList objects
        all the HDUList objects included in the validated input_data object

    """

    if isinstance(input_data, pf.HDUList):

        validated_fits_list = [input_data]

    elif isinstance(input_data, DataSet):

        validated_fits_list = input_data.spectroscopic_images

    else:
        raise PYWFC3InputError('Invalid input. Parameter input_data, should be either a HDUList or a DataSet object.')

    return validated_fits_list


def fits_list_size(input_data):

    if isinstance(input_data, pf.HDUList):
        return len(fits_list(input_data))
    elif isinstance(input_data, DataSet):
        if input_data.splitted:
            return len(fits_list(input_data)[0])
        else:
            return len(fits_list(input_data))


def save_figure(directory, figure=None, name=None, transparent=True):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    if not name:
        for fig in range(100000):
            if not os.path.isfile(directory + '/f{0}00.eps'.format(fig)):
                if not figure:
                    plt.savefig(directory + '/f{0}00.eps'.format(fig), bbox_inches='tight', transparent=transparent)
                    plt.close()
                else:
                    figure.savefig(directory + '/f{0}00.eps'.format(fig), bbox_inches='tight', transparent=transparent)
                break
    else:
        if not figure:
            plt.savefig(directory + '/{0}.eps'.format(name), bbox_inches='tight', transparent=transparent)
            plt.close()
        else:
            figure.savefig(directory + '/{0}.eps'.format(name), bbox_inches='tight', transparent=transparent)


def adjust_ticks():

    xlim1, xlim2 = plt.xlim()
    xticks = plt.xticks()[0]
    dxticks = xticks[1] - xticks[0]

    if xticks[1] - dxticks / 2 < xlim1:
        new_xlim1 = xticks[1] - dxticks / 2
        xticks = xticks[1:]
    else:
        new_xlim1 = xticks[0] - dxticks / 2

    if xticks[-2] + dxticks / 2 > xlim2:
        new_xlim2 = xticks[-2] + dxticks / 2
        xticks = xticks[:-1]
    else:
        new_xlim2 = xticks[-1] + dxticks / 2

    plt.xticks(xticks, xticks, fontsize=15)
    plt.xlim(new_xlim1, new_xlim2)

    ylim1, ylim2 = plt.ylim()
    yticks = plt.yticks()[0]
    dyticks = yticks[1] - yticks[0]

    if yticks[1] - dyticks / 2 < ylim1:
        new_ylim1 = yticks[1] - dyticks / 2
        yticks = yticks[1:]
    else:
        new_ylim1 = yticks[0] - dyticks / 2

    if yticks[-2] + dyticks / 2 > ylim2:
        new_ylim2 = yticks[-2] + dyticks / 2
        yticks = yticks[:-1]
    else:
        new_ylim2 = yticks[-1] + dyticks / 2

    plt.yticks(yticks, yticks, fontsize=15)
    plt.ylim(new_ylim1, new_ylim2)


# pipeline variables


class PipelineVariable():

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
                    if value == 'auto':
                        self.value = value
                    else:
                        try:
                            self.value = self.instance(value)
                        except ValueError:
                            raise PYWFC3InputError('Input {0} parameter is not valid, '
                                                   '{1} is expected.'.format(self.name, self.instance))
                else:
                    try:
                        self.value = self.instance(value)
                    except ValueError:
                        raise PYWFC3InputError('Input {0} parameter is not valid, '
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
                    raise PYWFC3InputError('Input {0} parameter is not valid, '
                                           '{1} is expected.'.format(self.name, self.instance))
            else:
                try:
                    self.value = self.instance(value)
                except ValueError:
                        raise PYWFC3InputError('Input {0} parameter is not valid, '
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

        if isinstance(dictionary, PipelineVariable):
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

        if isinstance(dictionary, PipelineVariable):
            if sub_dictionary is None:
                dictionary.value[self.keyword] = self.instance(self.value)
            else:
                dictionary.value[sub_dictionary][self.keyword] = self.instance(self.value)
        else:
            if sub_dictionary is None:
                dictionary[self.keyword] = self.instance(self.value)
            else:
                dictionary[sub_dictionary][self.keyword] = self.instance(self.value)

    def from_var(self, variable):

        if isinstance(variable, PipelineVariable):
            self.set(variable.value)
        else:
            self.set(variable)

    def custom(self, value=None):

        if value is None:
            return PipelineVariable(self.name, self.keyword, self.value, self.instance, self.kind,
                                    self.hdu_image, self.auto_fill)
        else:
            return PipelineVariable(self.name, self.keyword, value, self.instance, self.kind, self.hdu_image,
                                    self.auto_fill)

    def custom_from_fits(self, fits, position=0):

        x = self.custom()
        x.from_fits(fits, position)

        return x

    def custom_from_dictionary(self, dictionary, sub_dictionary=None):

        x = self.custom()
        x.from_dictionary(dictionary, sub_dictionary)

        return x


class PipelineVariablesSet():

    def __init__(self):

        # pipeline
        self.pipeline_files_location = PipelineVariable(
            'pipeline_files_location', value=os.path.abspath(os.path.dirname(__file__)), instance=str)
        self.pipeline_files_directory = PipelineVariable(
            'pipeline_files_directory', value='wfc3_calibration_files', instance=str)
        self.calibration_last_update = PipelineVariable(
            'calibration_last_update', value=20170410, instance=int)
        self.calibration_last_update_file = PipelineVariable(
            'calibration_last_update_file', value='last_update_calibration.txt', instance=str)
        self.calibration_zip_file = PipelineVariable(
            'calibration_zip_file', value='wfc3_calibration_files.zip', instance=str)
        self.calibration_directory = PipelineVariable(
            'calibration_directory', value='wfc3_calibration_files', instance=str)
        self.calibration_url = PipelineVariable(
            'calibration_url',
            value='http://zuserver2.star.ucl.ac.uk/~atsiaras/wfc3_calibration_files.zip', instance=str)

        # process visit
        self.data_directory = PipelineVariable(
            'data_directory', value='.', instance=str, kind='u')
        self.reduction = PipelineVariable(
            'reduction', value=True, instance=bool, kind='u')
        self.splitting = PipelineVariable(
            'splitting', value=False, instance=bool, kind='u')
        self.extraction = PipelineVariable(
            'extraction', value=True, instance=bool, kind='u')
        self.splitting_extraction = PipelineVariable(
            'splitting_extraction', value=False, instance=bool, kind='u')
        self.fitting_white = PipelineVariable(
            'fitting_white', value=True, instance=bool, kind='u')
        self.fitting_spectrum = PipelineVariable(
            'fitting_spectrum', value=True, instance=bool, kind='u')
        self.raw_data_directory = PipelineVariable(
            'raw_data_directory', value='raw_data', instance=str)
        self.reduced_data_directory = PipelineVariable(
            'reduced_data_directory', value='reduced_data', instance=str)
        self.splitted_data_directory = PipelineVariable(
            'splitted_data_directory', value='splitted_data', instance=str)
        self.figures_directory = PipelineVariable(
            'figures_directory', value='figures', instance=str)
        self.output_directory = PipelineVariable(
            'output_directory', value='results', instance=str)
        self.output_directory_copy = PipelineVariable(
            'output_directory_copy', value='results', instance=str, kind='u')
        self.light_curve_file = PipelineVariable(
            'light_curve_file', value='extracted_light_curves', instance=str)
        self.fitting_file = PipelineVariable(
            'fitting_file', value='fitting_results', instance=str)
        self.apply_bias = PipelineVariable(
            'apply_bias', value=True, instance=bool)
        self.apply_linearity = PipelineVariable(
            'apply_linearity', value=True, instance=bool)
        self.apply_dark = PipelineVariable(
            'apply_dark', value=True, instance=bool)
        self.apply_gain = PipelineVariable(
            'apply_gain', value=True, instance=bool)
        self.apply_sky = PipelineVariable(
            'apply_sky', value=True, instance=bool)
        self.apply_flat = PipelineVariable(
            'apply_flat', value=True, instance=bool)
        self.apply_bpcr = PipelineVariable(
            'apply_bpcr', value=True, instance=bool)

        # observation
        self.ra_target = PipelineVariable(
            'ra_target', keyword='RA_TARG', instance=float)
        self.dec_target = PipelineVariable(
            'dec_target', keyword='DEC_TARG', instance=float)
        self.observation_type = PipelineVariable(
            'observation_type', keyword='OBSTYPE', instance=str)
        self.grism = PipelineVariable(
            'grism', keyword='FILTER', instance=str)
        self.wfc3_aperture = PipelineVariable(
            'wfc3_aperture', keyword='APERTURE', instance=str)
        self.postarg1 = PipelineVariable(
            'postarg1', keyword='POSTARG1', instance=float)
        self.scan_direction = PipelineVariable(
            'scan_direction', keyword='POSTARG2', instance=float)
        self.exposure_start = PipelineVariable(
            'exposure_start', keyword='EXPSTART', instance=float)
        self.exposure_end = PipelineVariable(
            'exposure_end', keyword='EXPEND', instance=float)
        self.exposure_time = PipelineVariable(
            'exposure_time', keyword='EXPTIME', instance=float)
        self.sub_array_type = PipelineVariable(
            'sub_array_type', keyword='SUBTYPE', instance=str)
        self.sub_array_size = PipelineVariable(
            'sub_array_size', keyword='SUBARRAY', instance=int)
        self.total_samples = PipelineVariable(
            'total_samples', keyword='NSAMP', instance=int)
        self.sampling_sequence = PipelineVariable(
            'sampling_sequence', keyword='SAMP_SEQ', instance=str)
        self.sample_number = PipelineVariable(
            'sample_number', keyword='SAMPNUM', instance=int)
        self.sample_time = PipelineVariable(
            'sample_time', keyword='SAMPTIME', instance=float)

        # timing
        self.heliocentric_julian_date = PipelineVariable(
            'heliocentric_julian_date', keyword='HJD', instance=float)

        # bias and zero-read
        self.zero_read = PipelineVariable(
            'zero_read', keyword='ZREAD', instance=np.array, hdu_image=True)
        self.zero_read_error = PipelineVariable(
            'zero_read_error', keyword='ZREADE', instance=np.array, hdu_image=True)
        self.zero_read_flux = PipelineVariable(
            'zero_read_flux', keyword='ZFLUX', instance=np.array, hdu_image=True)
        self.zero_read_flux_error = PipelineVariable(
            'zero_read_flux_error', keyword='ZFLUXE', instance=np.array, hdu_image=True)
        self.reference_pixels_level = PipelineVariable(
            'reference_pixels_level', keyword='REFLEV', instance=float)

        # sky
        self.sky_detection_limit = PipelineVariable(
            'sky_detection_limit', keyword='SKYLIM', value=2.0, instance=float)
        self.sky_background_level = PipelineVariable(
            'sky_background_level', keyword='SKYBGLEV', instance=float)
        self.sky_area = PipelineVariable(
            'sky_area', keyword='SKYAREA', instance=np.array, hdu_image=True)

        # bpcr
        self.cr_neighbours = PipelineVariable(
            'cr_neighbours', keyword='CRXNB', value=6, instance=int)
        self.cr_detection_limit = PipelineVariable(
            'cr_detection_limit', keyword='CRLIM', value=3.0, instance=float)
        self.use_bpcr_fast_mode = PipelineVariable(
            'use_bpcr_fast_mode', keyword='FBPCR', value=True, instance=bool)
        self.bpcr_map = PipelineVariable(
            'bpcr_map', keyword='BPCRMAP', instance=np.array, hdu_image=True)

        # calibration
        self.comparison_index_forward = PipelineVariable(
            'comparison_index_forward', keyword='FCOMP', value=0, instance=int)
        self.comparison_index_reverse = PipelineVariable(
            'comparison_index_reverse', keyword='RCOMP', value=1, instance=int)
        self.target_x_offset = PipelineVariable(
            'target_x_offset', keyword='TARXOFF', value=0, instance=float, kind='u')
        self.target_y_offset = PipelineVariable(
            'target_y_offset', keyword='TARYOFF', value=0, instance=float, kind='u')
        self.use_standard_flat = PipelineVariable(
            'use_standard_flat', keyword='STDFLAT', value=True, instance=bool)

        self.spectrum_bottom = PipelineVariable(
            'spectrum_bottom', keyword='SPCBTTM', instance=int)
        self.spectrum_top = PipelineVariable(
            'spectrum_top', keyword='SPCTOP', instance=int)
        self.spectrum_left = PipelineVariable(
            'spectrum_left', keyword='SPCLEFT', instance=int)
        self.spectrum_right = PipelineVariable(
            'spectrum_right', keyword='SPCRIGHT', instance=int)
        self.spectrum_scan = PipelineVariable(
            'spectrum_scan', keyword='SPCSCAN', instance=bool)
        self.spectrum_direction = PipelineVariable(
            'spectrum_direction', keyword='SPECDIR', instance=int)
        self.first_spectrum_bottom = PipelineVariable(
            'first_spectrum_bottom', keyword='FSPCBTTM', instance=int)
        self.first_spectrum_top = PipelineVariable(
            'first_spectrum_top', keyword='FSPCTOP', instance=int)
        self.first_spectrum_scan = PipelineVariable(
            'first_spectrum_scan', keyword='FSPCSCAN', instance=bool)
        self.first_spectrum_direction = PipelineVariable(
            'first_spectrum_direction', keyword='FSPECDIR', instance=int)

        self.comparison_x_star = PipelineVariable(
            'comparison_x_star', keyword='CMPXSTAR', instance=float)
        self.x_star = PipelineVariable(
            'x_star', keyword='XSTAR', instance=float)
        self.x_shift = PipelineVariable(
            'x_shift', keyword='XSHIFT', instance=float)
        self.x_shift_error = PipelineVariable(
            'x_shift_error', keyword='XSHIFTE', instance=float)

        self.comparison_y_star = PipelineVariable(
            'comparison_y_star', keyword='CMPYSTAR', instance=float)
        self.y_star = PipelineVariable(
            'y_star', keyword='YSTAR', instance=float)
        self.y_shift = PipelineVariable(
            'y_shift', keyword='YSHIFT', instance=float)
        self.y_shift_error = PipelineVariable(
            'y_shift_error', keyword='YSHIFTE', instance=float)

        self.scan_length = PipelineVariable(
            'scan_length', keyword='LEN', instance=float)
        self.scan_length_error = PipelineVariable(
            'scan_length_error', keyword='LEN_ERR', instance=float)

        self.wdpt_constant_coefficient_1 = PipelineVariable(
            'wdpt_constant_coefficient_1', keyword='CSCOEFF1', instance=float)
        self.wdpt_constant_coefficient_2 = PipelineVariable(
            'wdpt_constant_coefficient_2', keyword='CSCOEFF2', instance=float)
        self.wdpt_constant_coefficient_3 = PipelineVariable(
            'wdpt_constant_coefficient_3', keyword='CSCOEFF3', instance=float)
        self.wdpt_slope_coefficient_1 = PipelineVariable(
            'wdpt_slope_coefficient_1', keyword='SLCOEFF1', instance=float)
        self.wdpt_slope_coefficient_2 = PipelineVariable(
            'wdpt_slope_coefficient_2', keyword='SLCOEFF2', instance=float)
        self.wdpt_slope_coefficient_3 = PipelineVariable(
            'wdpt_slope_coefficient_3', keyword='SLCOEFF3', instance=float)

        self.scan_frame = PipelineVariable(
            'scan_frame', keyword='SCANMAP', instance=np.array, hdu_image=True)
        self.wavelength_frame = PipelineVariable(
            'wavelength_frame', keyword='WMAP', instance=np.array, hdu_image=True)
        self.normalised_wavelength_frame = PipelineVariable(
            'normalised_wavelength_frame', keyword='NWMAP', instance=np.array, hdu_image=True)

        # extraction
        self.aperture_lower_extend = PipelineVariable(
            'aperture_lower_extend', value=-20.0, instance=float, kind='u')
        self.aperture_upper_extend = PipelineVariable(
            'aperture_upper_extend', value=20.0, instance=float, kind='u')
        self.extraction_method = PipelineVariable(
            'extraction_method', value='gauss', instance=str, kind='u')
        self.extraction_gauss_sigma = PipelineVariable(
            'extraction_gauss_sigma', value=45.7, instance=float, kind='u')
        self.white_lower_wavelength = PipelineVariable(
            'white_lower_wavelength', value=10800.0, instance=float, kind='u')
        self.white_upper_wavelength = PipelineVariable(
            'white_upper_wavelength', value=16800.0, instance=float, kind='u')
        self.bins_file = PipelineVariable(
            'bins_file', value='./bins_file.txt', instance=str, kind='u')
        self.bins_lower_wavelength = PipelineVariable('bins_lower_wavelength',
                                                      value=np.array([11240, 11440, 11640, 11840, 12040, 12240, 12440,
                                                                      12640, 12840, 13040, 13240, 13440, 13640, 13840,
                                                                      14040, 14240, 14440, 14640, 14840, 15040, 15240,
                                                                      15440, 15640, 15840, 16040, 16240]),
                                                      instance=np.array)
        self.bins_upper_wavelength = PipelineVariable('bins_upper_wavelength',
                                                      value=np.array([11440, 11640, 11840, 12040, 12240, 12440, 12640,
                                                                      12840, 13040, 13240, 13440, 13640, 13840, 14040,
                                                                      14240, 14440, 14640, 14840, 15040, 15240, 15440,
                                                                      15640, 15840, 16040, 16240, 16440]),
                                                      instance=np.array)
        self.heliocentric_julian_date_array = PipelineVariable(
            'heliocentric_julian_date_array', value=np.array([]), instance=np.array)
        self.spectrum_direction_array = PipelineVariable(
            'spectrum_direction_array', value=np.array([]), instance=np.array)
        self.sky_background_level_array = PipelineVariable(
            'sky_background_level_array', value=np.array([]), instance=np.array)
        self.x_star_array = PipelineVariable(
            'x_star_array', value=np.array([]), instance=np.array)
        self.x_shift_error_array = PipelineVariable(
            'x_shift_error_array', value=np.array([]), instance=np.array)
        self.y_star_array = PipelineVariable(
            'y_star_array', value=np.array([]), instance=np.array)
        self.y_shift_error_array = PipelineVariable(
            'y_shift_error_array', value=np.array([]), instance=np.array)
        self.scan_length_array = PipelineVariable(
            'scan_length_array', value=np.array([]), instance=np.array)
        self.scan_length_error_array = PipelineVariable(
            'scan_length_error_array', value=np.array([]), instance=np.array)
        self.lower_wavelength = PipelineVariable(
            'lower_wavelength', instance=float)
        self.upper_wavelength = PipelineVariable(
            'upper_wavelength', instance=float)
        self.flux_array = PipelineVariable(
            'flux', value=np.array([]), instance=np.array)
        self.error_array = PipelineVariable(
            'error', value=np.array([]), instance=np.array)
        self.ph_error_array = PipelineVariable(
            'ph_error', value=np.array([]), instance=np.array)
        self.white_dictionary = PipelineVariable(
            'white', value={}, instance=dict)
        self.bins_dictionaries = [PipelineVariable('bin_{0}'.format(str(ff + 1).zfill(2)), value={}, instance=dict)
                                  for ff in range(len(self.bins_lower_wavelength.value))]
        self.light_curve_split = PipelineVariable(
            'light_curve_split', 'split_', instance=str)

        # fitting
        self.apply_up_down_stream_correction = PipelineVariable('apply_up_down_stream_correction',
                                                                'USDS', True, bool, 'u')
        self.exclude_initial_orbits = PipelineVariable('exclude_initial_orbits', 'EXCL1', 1, int, 'u')
        self.exclude_final_orbits = PipelineVariable('exclude_final_orbits', 'EXCL2', 0, int, 'u')
        self.exclude_initial_orbit_points = PipelineVariable('exclude_initial_orbit_points', 'EXCP1', 0, int, 'u')
        self.white_ldc1 = PipelineVariable('white_ldc1', 'WHTLDC1', 0, float, 'u')
        self.white_ldc2 = PipelineVariable('white_ldc2', 'WHTLDC2', 0, float, 'u')
        self.white_ldc3 = PipelineVariable('white_ldc3', 'WHTLDC3', 0, float, 'u')
        self.white_ldc4 = PipelineVariable('white_ldc4', 'WHTLDC4', 0, float, 'u')
        self.fit_ldc1 = PipelineVariable('fit_ldc1', 'FITLDC1', False, bool, 'u')
        self.fit_ldc2 = PipelineVariable('fit_ldc2', 'FITLDC2', False, bool, 'u')
        self.fit_ldc3 = PipelineVariable('fit_ldc3', 'FITLDC3', False, bool, 'u')
        self.fit_ldc4 = PipelineVariable('fit_ldc4', 'FITLDC4', False, bool, 'u')
        self.bins_ldc1 = PipelineVariable('bins_ldc1', 'BINSLDC1', 'auto', np.array, auto_fill=True)
        self.bins_ldc2 = PipelineVariable('bins_ldc2', 'BINSLDC2', 'auto', np.array, auto_fill=True)
        self.bins_ldc3 = PipelineVariable('bins_ldc3', 'BINSLDC3', 'auto', np.array, auto_fill=True)
        self.bins_ldc4 = PipelineVariable('bins_ldc4', 'BINSLDC4', 'auto', np.array, auto_fill=True)
        self.planet = PipelineVariable('planet', 'PLANET', 'auto', str, 'u', auto_fill=True)
        self.method = PipelineVariable('method', 'METHOD', 'claret', str, 'u')
        self.rp_over_rs = PipelineVariable('rp_over_rs ', 'RP', 'auto', float, 'u', auto_fill=True)
        self.fp_over_fs = PipelineVariable('fp_over_fs', 'FP', 'auto', float, 'u', auto_fill=True)
        self.period = PipelineVariable('period', 'P', 'auto', float, 'u', auto_fill=True)
        self.sma_over_rs = PipelineVariable('sma_over_rs', 'A', 'auto', float, 'u', auto_fill=True)
        self.fit_sma_over_rs = PipelineVariable('fit_sma_over_rs', 'FITA', False, bool, 'u')
        self.eccentricity = PipelineVariable('eccentricity', 'E', 'auto', float, 'u', auto_fill=True)
        self.inclination = PipelineVariable('inclination', 'I', 'auto', float, 'u', auto_fill=True)
        self.fit_inclination = PipelineVariable('fit_inclination', 'FITI', False, bool, 'u')
        self.periastron = PipelineVariable('periastron', 'W', 'auto', float, 'u', auto_fill=True)
        self.mid_time = PipelineVariable('mid_time', 'MT', 'auto', float, 'u', auto_fill=True)
        self.fit_mid_time = PipelineVariable('fit_mid_time', 'FITMT', False, bool, 'u')
        self.second_order_ramp = PipelineVariable('second_order_ramp', 'RAMP2', False, bool, 'u')
        self.first_orbit_ramp = PipelineVariable('first_orbit_ramp', 'FOR', True, bool, 'u')
        self.mid_orbit_ramps = PipelineVariable('mid_orbit_ramps', 'MOR', True, bool, 'u')
        self.mcmc_iterations = PipelineVariable('mcmc_iterations', 'MCMCITER', 150000, int, 'u')
        self.mcmc_walkers = PipelineVariable('mcmc_walkers', 'MCMCWALK', 200, int, 'u')
        self.mcmc_burned_iterations = PipelineVariable('mcmc_burned_iterations', 'MCMCBURN', 50000, int, 'u')
        self.spectral_mcmc_iterations = PipelineVariable('spectral_mcmc_iterations', 'bla', 150000, int, 'u')
        self.spectral_mcmc_walkers = PipelineVariable('spectral_mcmc_walkers', 'blabla', 50, int, 'u')
        self.spectral_mcmc_burned_iterations = PipelineVariable(
            'spectral_mcmc_burned_iterations', 'blablabla', 50000, int, 'u')

    # def set_bins(self, bins_lower_wavelength, bins_upper_wavelength,
    #              bins_ldc1=None, bins_ldc2=None, bins_ldc3=None, bins_ldc4=None):
    #
    #     if bins_lower_wavelength is not None:
    #         if bins_upper_wavelength is not None:
    #
    #             self.bins_lower_wavelength.set(bins_lower_wavelength)
    #             self.bins_upper_wavelength.set(bins_upper_wavelength)
    #             self.bins_dictionaries = [PipelineVariable('bin_{0}'.format(str(ff + 1).zfill(2)),
    #                                                        value={}, instance=dict)
    #                                       for ff in range(len(bins_lower_wavelength))]
    #             if bins_ldc1 is None:
    #                 self.bins_ldc1.set('auto')
    #             else:
    #                 self.bins_ldc1.set(bins_ldc1)
    #             if bins_ldc2 is None:
    #                 self.bins_ldc2.set('auto')
    #             else:
    #                 self.bins_ldc2.set(bins_ldc2)
    #             if bins_ldc3 is None:
    #                 self.bins_ldc3.set('auto')
    #             else:
    #                 self.bins_ldc3.set(bins_ldc3)
    #             if bins_ldc4 is None:
    #                 self.bins_ldc4.set('auto')
    #             else:
    #                 self.bins_ldc4.set(bins_ldc4)

    def from_parameters_file(self, parameters_file):

        parameters_file = os.path.abspath(parameters_file)
        if not os.path.isfile(parameters_file):
            raise PYWFC3FileError('No such file: ' + parameters_file)
        print 'Loading parameters file: {} ...'.format(os.path.split(parameters_file)[1])

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
                        print 'WARNING: Parameter not in file, setting {0}={1}'.format(i, vars(self)[i].default_value)

        bins_file = os.path.abspath(self.bins_file.value)
        if os.path.isfile(bins_file):
            print 'Loading bins file: {} ...'.format(os.path.split(bins_file)[1])

            if len(np.loadtxt(bins_file, unpack=True)) == 2:
                bins_lower_wavelength, bins_upper_wavelength = np.loadtxt(bins_file, unpack=True)
                bins_ldc1 = None
                bins_ldc2 = None
                bins_ldc3 = None
                bins_ldc4 = None

            else:
                bins_lower_wavelength, bins_upper_wavelength, bins_ldc1, bins_ldc2, bins_ldc3, bins_ldc4 = \
                    np.loadtxt(bins_file, unpack=True)

            self.bins_lower_wavelength.set(bins_lower_wavelength)
            self.bins_upper_wavelength.set(bins_upper_wavelength)
            self.bins_dictionaries = [PipelineVariable('bin_{0}'.format(str(ff + 1).zfill(2)), value={}, instance=dict)
                                      for ff in range(len(bins_lower_wavelength))]
            if bins_ldc1 is None:
                self.bins_ldc1.set('auto')
            else:
                self.bins_ldc1.set(bins_ldc1)
            if bins_ldc2 is None:
                self.bins_ldc2.set('auto')
            else:
                self.bins_ldc2.set(bins_ldc2)
            if bins_ldc3 is None:
                self.bins_ldc3.set('auto')
            else:
                self.bins_ldc3.set(bins_ldc3)
            if bins_ldc4 is None:
                self.bins_ldc4.set('auto')
            else:
                self.bins_ldc4.set(bins_ldc4)

        else:
            print 'WARNING: No bins file found, defaults bins will be used.'

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


pipeline_variables = PipelineVariablesSet()


# calibration variables


class CalibrationVariable():

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

        if isinstance(input_data, DataSet):
            if input_data.splitted:
                input_data = fits_list(input_data)[0][0]
            else:
                input_data = fits_list(input_data)[0]

        calibration_directory_path = os.path.join(pipeline_variables.pipeline_files_location.value,
                                                  pipeline_variables.calibration_directory.value)

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
                    sample_numbers = [input_data[pp].header[pipeline_variables.sample_number.keyword]
                                      for pp in sci(input_data)]

                    for i in sci(selected_calibration_data):
                        if (selected_calibration_data[i].header[
                                pipeline_variables.sample_number.keyword] in sample_numbers
                            or selected_calibration_data[i].header[
                                pipeline_variables.sample_number.keyword] == min(sample_numbers) - 1):

                            dark_list = np.append(dark_list, selected_calibration_data[i])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 1])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 2])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 3])
                            dark_list = np.append(dark_list, selected_calibration_data[i + 4])

                    return pf.HDUList(list(dark_list))

                else:

                    array = (np.ones_like(selected_calibration_data[self.keyword].data) *
                             selected_calibration_data[self.keyword].data)

                    crop1 = len(np.array(array)) / 2 - len(input_data[sci(input_data)[0]].data) / 2
                    crop2 = len(np.array(array)) / 2 + len(input_data[sci(input_data)[0]].data) / 2

                    return array[crop1:crop2, crop1:crop2]

            else:
                return (np.ones_like(selected_calibration_data[0].header[self.keyword]) *
                        selected_calibration_data[0].header[self.keyword])

        else:
            raise PYWFC3LibraryError('No suitable dark current calibration file found, \n'
                                     'you must update your library manually. Search at: \n\n'
                                     '\thttps://hst-crds.stsci.edu/ \n\n'
                                     'for the latest dark current calibration file for\n\n'
                                     '\tSAMP_SEQ: {0} \n\tSUBTYTE: {1}'.format(input_data[0].header['SAMP_SEQ'],
                                                                               input_data[0].header['SUBTYPE']))


class CalibrationVariablesSet():

    def __init__(self):

        files_location = pipeline_variables.pipeline_files_location.value
        files_directory_path = os.path.join(files_location, pipeline_variables.pipeline_files_directory.value)
        calibration_last_update_file_path = os.path.join(files_directory_path,
                                                         pipeline_variables.calibration_last_update_file.value)
        calibration_zip_file_path = os.path.join(files_location, pipeline_variables.calibration_zip_file.value)
        calibration_directory_path = files_directory_path

        # update calibration files

        calibration_update = False
        if not os.path.isdir(files_directory_path):
            calibration_update = True
        elif not os.path.isfile(calibration_last_update_file_path):
            calibration_update = True
        elif (int(open(calibration_last_update_file_path).readlines()[0]) <
                pipeline_variables.calibration_last_update.value):
            calibration_update = True

        if calibration_update:
            try:
                print '\nDownloading calibration files...'

                if os.path.isdir(calibration_directory_path):
                    os.system('rm -rf {0}'.format(calibration_directory_path))

                urllib.urlretrieve(pipeline_variables.calibration_url.value, calibration_zip_file_path)

                os.system('unzip {0} -d {1}{2}'.format(calibration_zip_file_path, files_location, os.sep))
                os.system('rm {0}'.format(calibration_zip_file_path))
                os.system('rm -rf {0}{1}__MACOSX'.format(files_location, os.sep))

            except IOError:
                raise PYWFC3LibraryError('Failed to update wfc3 calibration files.')

        self.super_zero_read = CalibrationVariable(
            'lin_custom', 'ZSCI', hdu_image=True)
        self.super_zero_read_error = CalibrationVariable(
            'lin_custom', 'ZERR', hdu_image=True)
        self.ccd_gain = CalibrationVariable(
            'ccd_custom', 'CCDGAIN', hdu_image=True)
        self.ccd_read_noise = CalibrationVariable(
            'ccd_custom', 'CCDREADN', hdu_image=True)
        self.linearity_coefficient_1 = CalibrationVariable(
            'lin_custom', 'LINC1', hdu_image=True)
        self.linearity_coefficient_2 = CalibrationVariable(
            'lin_custom', 'LINC2', hdu_image=True)
        self.linearity_coefficient_3 = CalibrationVariable(
            'lin_custom', 'LINC3', hdu_image=True)
        self.linearity_coefficient_4 = CalibrationVariable(
            'lin_custom', 'LINC4', hdu_image=True)
        self.super_dark = CalibrationVariable(
            'drk', 'SCI', dark_mode=True, hdu_image=True,
            matching_keywords=[pipeline_variables.sampling_sequence.keyword, pipeline_variables.sub_array_type.keyword])
        self.mean_gain = CalibrationVariable(
            'ccd_custom', 'MGAIN')
        self.gain_profile = CalibrationVariable(
            'pfl', 'SCI', hdu_image=True)
        self.master_sky = CalibrationVariable(
            'sky', 'PRIMARY', hdu_image=True, matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_coefficient_1 = CalibrationVariable(
            'flat', 'FFC1', hdu_image=True, matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_coefficient_2 = CalibrationVariable(
            'flat', 'FFC2', hdu_image=True, matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_coefficient_3 = CalibrationVariable(
            'flat', 'FFC3', hdu_image=True, matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_coefficient_4 = CalibrationVariable(
            'flat', 'FFC4', hdu_image=True, matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_min_wavelength = CalibrationVariable(
            'flat', 'WMIN', matching_keywords=[pipeline_variables.grism.keyword])
        self.flat_field_max_wavelength = CalibrationVariable(
            'flat', 'WMAX', matching_keywords=[pipeline_variables.grism.keyword])
        self.bad_pixels = CalibrationVariable(
            'bpx_custom', 'PRIMARY', hdu_image=True)
        self.trace_at0 = CalibrationVariable(
            'calib', 'AT0', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_at1 = CalibrationVariable(
            'calib', 'AT1', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_at2 = CalibrationVariable(
            'calib', 'AT2', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_at3 = CalibrationVariable(
            'calib', 'AT3', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_at4 = CalibrationVariable(
            'calib', 'AT4', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_at5 = CalibrationVariable(
            'calib', 'AT5', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_bt0 = CalibrationVariable(
            'calib', 'BT0', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_bt1 = CalibrationVariable(
            'calib', 'BT1', matching_keywords=[pipeline_variables.grism.keyword])
        self.trace_bt2 = CalibrationVariable(
            'calib', 'BT2', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw0 = CalibrationVariable(
            'calib', 'AW0', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw1 = CalibrationVariable(
            'calib', 'AW1', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw2 = CalibrationVariable(
            'calib', 'AW2', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw3 = CalibrationVariable(
            'calib', 'AW3', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw4 = CalibrationVariable(
            'calib', 'AW4', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_aw5 = CalibrationVariable(
            'calib', 'AW5', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_bw0 = CalibrationVariable(
            'calib', 'BW0', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_bw1 = CalibrationVariable(
            'calib', 'BW1', matching_keywords=[pipeline_variables.grism.keyword])
        self.wsol_bw2 = CalibrationVariable(
            'calib', 'BW2', matching_keywords=[pipeline_variables.grism.keyword])


calibration_variables = CalibrationVariablesSet()


# data set object


class DataSet():

    def __init__(self, data_set_directory_path=None):

        if data_set_directory_path is None:

            self.spectroscopic_images = []
            self.splitted = False

        else:

            sub_directories = glob.glob(
                data_set_directory_path + '/' + pipeline_variables.reduced_data_directory.value + '*/')

            if len(sub_directories) == 0:

                nsamp = []
                final_list = []
                direct_image = False

                files = glob.glob(os.path.join(data_set_directory_path, '*.fits'))
                for i in files:
                    j = pf.open(i, memmap=False)
                    if j[0].header[pipeline_variables.observation_type.keyword] == 'SPECTROSCOPIC':
                        final_list.append([j[0].header[pipeline_variables.exposure_start.keyword],
                                           os.path.split(i)[1], j])
                        nsamp.append(j[0].header[pipeline_variables.total_samples.keyword])
                    elif not direct_image:
                        direct_image = pf.open(i, memmap=False)

                nsamps = [int(np.median(np.array(nsamp)))]

                final_list.sort()
                list_of_times, list_of_files, list_of_fits = np.swapaxes(final_list, 0, 1)

                outliers = True
                while outliers:
                    outliers = False
                    for i in range(len(list_of_fits)):
                        if list_of_fits[i][0].header[pipeline_variables.total_samples.keyword] not in nsamps:
                            list_of_fits = np.delete(list_of_fits, i)
                            list_of_files = np.delete(list_of_files, i)
                            outliers = True
                            break

                if not direct_image:
                    raise PYWFC3InputError('Direct image not included in data set.')

                if len(list_of_fits) < 2:
                    raise PYWFC3InputError('A data set should contain more than one spectroscopic images.')

                self.file_names = list_of_files
                self.spectroscopic_images = list_of_fits
                self.direct_image = direct_image
                self.splitted = False

            else:
                self.file_names = []
                self.spectroscopic_images = []
                self.direct_image = []
                self.splitted = True
                self._data_set_directory_path = data_set_directory_path

                for data_set_directory_path in sub_directories:

                    nsamp = []
                    final_list = []
                    direct_image = False

                    files = glob.glob(os.path.join(data_set_directory_path, '*.fits'))
                    for i in files:
                        j = pf.open(i, memmap=False)
                        if j[0].header[pipeline_variables.observation_type.keyword] == 'SPECTROSCOPIC':
                            final_list.append([j[0].header[pipeline_variables.exposure_start.keyword],
                                               os.path.split(i)[1], j])
                            nsamp.append(j[0].header[pipeline_variables.total_samples.keyword])
                        elif not direct_image:
                            direct_image = pf.open(i, memmap=False, mode='update')

                    nsamps = [int(np.median(np.array(nsamp)))]

                    final_list.sort()
                    list_of_times, list_of_files, list_of_fits = np.swapaxes(final_list, 0, 1)

                    outliers = True
                    while outliers:
                        outliers = False
                        for i in range(len(list_of_fits)):
                            if list_of_fits[i][0].header[pipeline_variables.total_samples.keyword] not in nsamps:
                                list_of_fits = np.delete(list_of_fits, i)
                                list_of_files = np.delete(list_of_files, i)
                                outliers = True
                                break

                    if not direct_image:
                        raise PYWFC3InputError('Direct image not included in data set.')

                    if len(list_of_fits) < 2:
                        raise PYWFC3InputError('A data set should contain more than one spectroscopic images.')

                    self.file_names = list_of_files
                    self.spectroscopic_images.append(list_of_fits)
                    self.direct_image = direct_image

    def save(self, export_directory, arrange=True, export_pipeline_variables_file='pipeline_variables.txt'):

        if os.path.isdir(export_directory):
            backup = export_directory + '_' + time.strftime('%y-%m-%d_%H-%M-%S')
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
            tools.fits_like(self.spectroscopic_images[i]).writeto(os.path.join(export_directory,
                                                                               self.file_names[i]))

        tools.fits_like(self.direct_image).writeto(os.path.join(export_directory, 'direct_image.fits'))

        if export_pipeline_variables_file:

            pipeline_variables.save(os.path.join(export_directory, export_pipeline_variables_file))

    def copy_split(self, split_number):

        x = DataSet()
        x.spectroscopic_images = self.spectroscopic_images[split_number]

        return x


# pipeline counter


class PipelineCounter():

    def __init__(self, task, total_iterations, show_every=1):

        self.task = task + '.' * (15 - len(task))
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

        out_of = ' ' * (len(str(self.total_iterations)) - len(str(self.current_iteration))) + \
                 str(self.current_iteration) + ' / ' + str(self.total_iterations)

        delta_time = time.time() - self.start_time

        time_left = str(datetime.timedelta(
            seconds=int((self.total_iterations - self.current_iteration) * delta_time / self.current_iteration)))

        total_time = str(datetime.timedelta(seconds=int(delta_time)))
            
        if int(self.show):

            sys.stdout.write('\r\033[K')
            sys.stdout.write(self.task + ': ' + out_of + '   time left: ' + time_left + '   total time: ' + total_time)
            sys.stdout.flush()
            self.show = 0

        if self.current_iteration == self.total_iterations and self.total_iterations > 1:
            print ''
