
import os
import glob
import time
import shutil
import pylightcurve as plc

from astropy.io import fits
from urllib.request import urlretrieve

from iraclis import __version__
from iraclis.errors import *
from iraclis.classes import *
from iraclis.images1_reduction import *
from iraclis.images2_calibration import *
from iraclis.images3_photometry import *
from iraclis.lightcurves1_fitting import *


def process_visit(parameters_file=None, data_directory=None, procedure=None, par_string=None):

    print('')
    print('Iraclis {0}'.format(__version__))

    # user inputs

    variables.from_parameters_file(parameters_file, data_directory)

    # for developers use only

    variables.overwrite(procedure, par_string)

    # load the parameters file

    directory = os.path.abspath(variables.data_directory.value)

    if not os.path.isdir(directory):
        raise IraclisFileError('No such directory: ' + directory)
    else:
        print('')
        print('Processing directory: ', os.path.split(directory)[1])
        print('                   @: ', os.path.split(directory)[0])
        print('')

    # set paths

    raw_data_directory = os.path.join(directory, variables.raw_data_directory.value)
    reduced_data_directory = os.path.join(directory, variables.reduced_data_directory.value)
    output_directory = os.path.join(directory, variables.output_directory.value)
    output_directory_copy = os.path.join(directory, variables.output_directory_copy.value)
    figures_directory = variables.figures_directory.value
    extraction_figures_directory = os.path.join(output_directory, '{0}_extraction'.format(figures_directory))
    fitting_figures_directory = os.path.join(output_directory, '{0}_fitting'.format(figures_directory))
    lc_file = os.path.join(output_directory, '{0}.pickle'.format(variables.light_curve_file.value))
    fit_file = os.path.join(output_directory, '{0}.pickle'.format(variables.fitting_file.value))
    splitted_data_directory = os.path.join(directory, variables.splitted_data_directory.value)
    splitted_reduced_data_directory = os.path.join(splitted_data_directory,
                                                   variables.reduced_data_directory.value)

    # set process steps and detect inconsistencies

    if variables.fitting_spectrum.value and not variables.fitting_white.value:
        if not os.path.isfile(fit_file):
            print('\nResetting Fitting-white to True ...')
            variables.fitting_white.set(True)

    if (variables.fitting_white.value and not variables.extraction.value
       and not variables.splitting_extraction.value):
        if not os.path.isfile(lc_file):
            print('\nResetting Extraction to True ...')
            variables.extraction.set(True)

    if variables.extraction.value and not variables.reduction.value:
        if not os.path.isdir(reduced_data_directory):
            print('\nResetting Reduction to True ...')
            variables.reduction.set(True)

    if variables.splitting_extraction.value and not variables.splitting.value:
        if not os.path.isdir(splitted_data_directory):
            print('\nResetting Splitting to True ...')
            variables.splitting.set(True)

    if variables.splitting_extraction.value and variables.extraction.value:
            print('\nResetting Splitting extraction to False ...')
            variables.splitting_extraction.set(False)

    if variables.splitting.value and not variables.reduction.value:
        if not os.path.isdir(reduced_data_directory):
            print('\nResetting Reduction to True ...')
            variables.reduction.set(True)

    if variables.reduction.value:
        if not os.path.isdir(raw_data_directory):
            print('\nMoving raw images in {} ...'.format(variables.raw_data_directory.value))
            os.mkdir(raw_data_directory)
            for fits_file in glob.glob(os.path.join(directory, '*.fits')):
                shutil.move(fits_file, os.path.join(raw_data_directory, os.path.split(fits_file)[1]))

    # reduction and calibration

    if variables.reduction.value:

        # load raw images

        print('Loading raw images ...')
        data_set = DataSet(raw_data_directory)

        data_set = timing(data_set)
        data_set = bias(data_set)
        data_set = linearity(data_set)
        data_set = dark(data_set)
        data_set = gain(data_set)
        data_set = sky(data_set)
        data_set = calibration(data_set)
        data_set = flat(data_set)
        data_set = bpcr(data_set)

        # save reduced images

        print('Saving reduced frames in {} ...'.format(os.path.split(reduced_data_directory)[1]))
        data_set.save(reduced_data_directory)
        del data_set
        variables.save(os.path.join(reduced_data_directory, 'parameters.txt'))
        print('')

    # splitting

    if variables.splitting.value:

        print('Loading raw images ...')
        raw_data_set = DataSet(raw_data_directory)

        # load reduced images

        print('Loading reduced images ...')
        data_set = DataSet(reduced_data_directory)

        # split images

        if os.path.isdir(splitted_data_directory):
            backup = '{0}_{1}'.format(splitted_data_directory, time.strftime('%y-%m-%d_%H-%M-%S'))
            shutil.copytree(splitted_data_directory, backup)
            shutil.rmtree(splitted_data_directory)

        os.mkdir(splitted_data_directory)

        total_samples = raw_data_set.spectroscopic_images[0][0].header[variables.total_samples.keyword] - 1

        for sample in range(total_samples):

            print('Splitting sample {0}:'.format(sample + 1))

            for i in range(len(raw_data_set.spectroscopic_images)):

                frame = plc.copy_fits(data_set.spectroscopic_images[i])
                frame2 = plc.copy_fits(raw_data_set.spectroscopic_images[i])
                new_frame = [fits.PrimaryHDU(header=frame[0].header, data=frame[0].data),
                             fits.ImageHDU(header=frame2[1 + sample * 5].header, data=frame2[1 + sample * 5].data),
                             fits.ImageHDU(header=frame2[2 + sample * 5].header, data=frame2[2 + sample * 5].data),
                             fits.ImageHDU(header=frame2[3 + sample * 5].header, data=frame2[3 + sample * 5].data),
                             fits.ImageHDU(header=frame2[4 + sample * 5].header, data=frame2[4 + sample * 5].data),
                             fits.ImageHDU(header=frame2[5 + sample * 5].header, data=frame2[5 + sample * 5].data),
                             fits.ImageHDU(header=frame2[6 + sample * 5].header, data=frame2[6 + sample * 5].data),
                             fits.ImageHDU(header=frame2[7 + sample * 5].header, data=frame2[7 + sample * 5].data),
                             fits.ImageHDU(header=frame2[8 + sample * 5].header, data=frame2[8 + sample * 5].data),
                             fits.ImageHDU(header=frame2[9 + sample * 5].header, data=frame2[9 + sample * 5].data),
                             fits.ImageHDU(header=frame2[10 + sample * 5].header, data=frame2[10 + sample * 5].data),
                             fits.ImageHDU(header=frame['SKYAREA'].header, data=frame['SKYAREA'].data,
                                           name='SKYAREA'),
                             fits.ImageHDU(header=frame['SCANMAP'].header, data=frame['SCANMAP'].data,
                                           name='SCANMAP'),
                             fits.ImageHDU(header=frame['WMAP'].header, data=frame['WMAP'].data, name='WMAP'),
                             fits.ImageHDU(header=frame['NWMAP'].header, data=frame['NWMAP'].data, name='NWMAP'),
                             fits.ImageHDU(header=frame['BPCRMAP'].header, data=frame['BPCRMAP'].data,
                                           name='BPCRMAP'),
                             ]

                data_set.spectroscopic_images[i] = fits.HDUList(new_frame)

            data_set = bias(data_set)
            data_set = linearity(data_set)
            data_set = dark(data_set, splitting=True)
            data_set = gain(data_set)
            data_set = sky(data_set, splitting=True)
            data_set = calibration(data_set, splitting=True)
            data_set = flat(data_set)
            data_set = bpcr(data_set, splitting=True)
            data_set.save('{0}_{1}'.format(splitted_reduced_data_directory, str(sample + 1).zfill(2)))

        del data_set
        del raw_data_set

        print('')

    # create a backup if output directory exists

    if (variables.extraction.value or variables.splitting_extraction.value or
            variables.fitting_white.value or variables.fitting_spectrum.value):

        if os.path.isdir(output_directory):
            backup = '{0}_{1}'.format(output_directory, time.strftime('%y-%m-%d_%H-%M-%S'))
            shutil.copytree(output_directory, backup)
        else:
            os.mkdir(output_directory)

    # extraction of light curves

    if variables.extraction.value or variables.splitting_extraction.value:

        if variables.extraction.value:

            # load reduced images

            print('Loading reduced images ...')
            data_set = DataSet(reduced_data_directory)

            # use the photometry function to extract the light curves

            data_set, light_curve = photometry(data_set)

        else:

            print('Loading reduced images ...')
            data_set = DataSet(reduced_data_directory)

            # use the photometry function to extract the light curves

            data_set, light_curve = photometry(data_set)

            star_y_position_array = variables.y_star_array.custom_from_dictionary(light_curve)
            y_shift_error_array = variables.y_shift_error_array.custom_from_dictionary(light_curve)
            scan_length_array = variables.scan_length_array.custom_from_dictionary(light_curve)
            scan_length_error_array = variables.scan_length_error_array.custom_from_dictionary(light_curve)

            # load reduced splitted images

            print('Loading splitted reduced images ...')
            data_set = DataSet(splitted_data_directory)

            # use the split_photometry function to extract the light curves

            data_set, light_curve = split_photometry(data_set)

            star_y_position_array.to_dictionary(light_curve)
            y_shift_error_array.to_dictionary(light_curve)
            scan_length_array.to_dictionary(light_curve)
            scan_length_error_array.to_dictionary(light_curve)

        # save extraction results

        print('Saving extracted light-curves in {} ...'.format(str(os.sep).join(lc_file.split(os.sep)[-2:])))
        plc.save_dict(light_curve, lc_file)

        # plot extraction diagnostics and results

        print('Saving extraction plots in {} ...'.format(os.path.split(extraction_figures_directory)[1]))
        plot_photometry(data_set, light_curve, extraction_figures_directory)

        del data_set

        variables.save(os.path.join(output_directory, 'parameters.txt'))

        print('')

    # fitting the white and the spectral light curves

    if variables.fitting_white.value:

        light_curve = plc.open_dict(lc_file)
        light_curve = fitting(light_curve, fitting_spectrum=False)

        print('Saving fitting results in {} ...'.format(str(os.sep).join(fit_file.split(os.sep)[-2:])))
        plc.save_dict(light_curve, fit_file)

        print('Saving fitting plots in {} ...'.format(os.path.split(fitting_figures_directory)[1]))
        plot_fitting(light_curve, fitting_figures_directory)

        variables.save(os.path.join(output_directory, 'parameters.txt'))

        print('')

    # fitting the spectral light curves

    if variables.fitting_spectrum.value:

        light_curve = plc.open_dict(lc_file)
        fitted_white_light_curve = plc.open_dict(fit_file)
        light_curve = fitting(light_curve, fitted_white_light_curve=fitted_white_light_curve,
                              fitting_spectrum=variables.fitting_spectrum.value)

        print('Saving fitting results in {} ...'.format(str(os.sep).join(fit_file.split(os.sep)[-2:])))
        plc.save_dict(light_curve, fit_file)

        print('Saving fitting plots in {} ...'.format(os.path.split(fitting_figures_directory)[1]))
        plot_fitting(light_curve, fitting_figures_directory)

        variables.save(os.path.join(output_directory, 'parameters.txt'))

        print('')

    # copy results
    # create a backup if output directory exists

    if variables.output_directory_copy.value != 'False':

        if output_directory_copy != output_directory:

            if os.path.isdir(output_directory_copy):
                backup = '{0}_{1}'.format(output_directory_copy, time.strftime('%y-%m-%d_%H-%M-%S'))
                shutil.copytree(output_directory_copy, backup)
                shutil.rmtree(output_directory_copy)

            print('Copying {} to {} ...'.format(os.path.split(output_directory)[1],
                                                os.path.split(output_directory_copy)[1]))
            if os.path.isdir(output_directory):
                shutil.copytree(output_directory, output_directory_copy)

            print('')

    print('Processed directory: ', os.path.split(directory)[1])
    print('                  @: ', os.path.split(directory)[0])
    print('')


def run_test():
    dataset_files = [
        'icy021ljq_raw.fits',
        'icy021l6q_raw.fits',
        'icy021m2q_raw.fits',
        'icy021kuq_raw.fits',
        'icy021lgq_raw.fits',
        'icy021kmq_raw.fits',
        'icy021k1q_raw.fits',
        'icy021kxq_raw.fits',
        'icy021lrq_raw.fits',
        'icy021l0q_raw.fits',
        'icy021lyq_raw.fits',
        'icy021ksq_raw.fits',
        'icy021kfq_raw.fits',
        'icy021llq_raw.fits',
        'icy021jzq_raw.fits',
        'icy021ltq_raw.fits',
        'icy021kkq_raw.fits',
        'icy021laq_raw.fits',
        'icy021lkq_raw.fits',
        'icy021kaq_raw.fits',
        'icy021m3q_raw.fits',
        'icy021ktq_raw.fits',
        'icy021l7q_raw.fits',
        'icy021lfq_raw.fits',
        'icy021klq_raw.fits',
        'icy021kyq_raw.fits',
        'icy021lsq_raw.fits',
        'icy021k0q_raw.fits',
        'icy021lxq_raw.fits',
        'icy021kgq_raw.fits',
        'icy021lmq_raw.fits',
        'icy021luq_raw.fits',
        'icy021k6q_raw.fits',
        'icy021kjq_raw.fits',
        'icy021l4q_raw.fits',
        'icy021m0q_raw.fits',
        'icy021kwq_raw.fits',
        'icy021lhq_raw.fits',
        'icy021kbq_raw.fits',
        'icy021k3q_raw.fits',
        'icy021l9q_raw.fits',
        'icy021kzq_raw.fits',
        'icy021lpq_raw.fits',
        'icy021leq_raw.fits',
        'icy021koq_raw.fits',
        'icy021kdq_raw.fits',
        'icy021lnq_raw.fits',
        'icy021l2q_raw.fits',
        'icy021k8q_raw.fits',
        'icy021kqq_raw.fits',
        'icy021kiq_raw.fits',
        'icy021lcq_raw.fits',
        'icy021k5q_raw.fits',
        'icy021lvq_raw.fits',
        'icy021m1q_raw.fits',
        'icy021kvq_raw.fits',
        'icy021l5q_raw.fits',
        'icy021liq_raw.fits',
        'icy021kcq_raw.fits',
        'icy021lqq_raw.fits',
        'icy021k2q_raw.fits',
        'icy021l8q_raw.fits',
        'icy021ldq_raw.fits',
        'icy021knq_raw.fits',
        'icy021keq_raw.fits',
        'icy021loq_raw.fits',
        'icy021jxq_flt.fits',
        'icy021lzq_raw.fits',
        'icy021kpq_raw.fits',
        'icy021l3q_raw.fits',
        'icy021k9q_raw.fits',
        'icy021lbq_raw.fits',
        'icy021lwq_raw.fits',
        'icy021jyq_raw.fits',
        'icy021k4q_raw.fits',
    ]

    desktop = os.path.join(os.path.expanduser('~'), 'Desktop')
    if not os.path.isdir(desktop):
        desktop = glob.glob(os.path.join(os.path.expanduser('~'), '*', 'Desktop'))[0]
    if not os.path.isdir(desktop):
        desktop = os.path.expanduser('~')

    destination = os.path.join(desktop, 'iraclis_test_dataset_hatp26b')

    if not os.path.isdir(destination):
        os.mkdir(destination)

    if not os.path.isdir(os.path.join(destination, 'raw_data')):
        for num, dataset_file in enumerate(dataset_files):
            print('{0}/{1}: '.format(num + 1, len(dataset_files)), dataset_file)
            if not os.path.isfile(os.path.join(destination, dataset_file)):
                urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(
                    dataset_file), os.path.join(destination, dataset_file))

    parameters_file = os.path.join(destination, 'iraclis_test_dataset_hatp26b_parameters.txt')
    w = open(parameters_file, 'w')
    w.write('\n'.join([
        'data_directory                     {0}          '.format(destination),
        '# directory path                                ',
        'output_directory_copy              False        ',
        '# directory name/False                          ',
        '                                                ',
        'reduction                          True         ',
        '# True/False                                    ',
        'splitting                          False        ',
        '# True/False                                    ',
        'extraction                         True         ',
        '# True/False                                    ',
        'splitting_extraction               False        ',
        '# True/False                                    ',
        'fitting_white                      True         ',
        '# True/False                                    ',
        'fitting_spectrum                   True         ',
        '# True/False                                    ',
        '                                                ',
        'target_x_offset                    0.0          ',
        '# number                                        ',
        'target_y_offset                    0.0          ',
        '# number                                        ',
        'aperture_lower_extend              -20.0        ',
        '# number                                        ',
        'aperture_upper_extend              20.0         ',
        '# number                                        ',
        'extraction_method                  gauss        ',
        '# gauss/integral                                ',
        'extraction_gauss_sigma             47.0         ',
        '# number                                        ',
        '                                                ',
        'method                             claret       ',
        '# claret/linear/quad/sqrt                       ',
        '                                                ',
        'white_lower_wavelength             default      ',
        '# number/default                                ',
        'white_upper_wavelength             default      ',
        '# number/default                                ',
        'white_ldc1                         default      ',
        '# number/default                                ',
        'white_ldc2                         default      ',
        '# number/default                                ',
        'white_ldc3                         default      ',
        '# number/default                                ',
        'white_ldc4                         default      ',
        '# number/default                                ',
        '                                                ',
        '# Comment: When the above parameters are set '
        'to default, the white light-curve limits will '
        'be 10880.0 - 16800.0 Angstroms for G141 and '
        '8000.0 - 11500.0 Angstroms for G102.            ',
        '                                                ',
        'bins_file                          default_low  ',
        '# file path/default_high/default_low/default_'
        'vlow                                            ',
        '                                                ',
        '# Comment: A bins file should contain 4 '
        'columns: '
        '1 - bins lower edge in Angstrom, '
        '2 - bins upper edge in Angstrom, '
        '3 - first limb darkening coefficient, '
        '4 - second limb darkening coefficient, '
        '5 - third limb darkening coefficient, '
        '6 - forth limb darkening coefficient. '
        'An example of a bins file can be found in '
        'iraclis_test_dataset_hatp26b_bins.txt           ',
        '                                                ',
        'planet                             auto         ',
        '# name - no spaces/auto                         ',
        'star_teff                          auto         ',
        '# number/auto                                   ',
        'star_logg                          auto         ',
        '# number/auto                                   ',
        'star_meta                          auto         ',
        '# number/auto                                   ',
        'rp_over_rs                         auto         ',
        '# number/auto                                   ',
        'fp_over_fs                         auto         ',
        '# number/auto                                   ',
        'period                             auto         ',
        '# number/auto                                   ',
        'sma_over_rs                        auto         ',
        '# number/auto                                   ',
        'eccentricity                       auto         ',
        '# number/auto                                   ',
        'inclination                        auto         ',
        '# number/auto                                   ',
        'periastron                         auto         ',
        '# number/auto                                   ',
        'mid_time                           auto         ',
        '# number/auto                                   ',
        '                                                ',
        '# Comment: You can set any of the above 12 '
        'parameters to auto, to use the data from '
        'a built-in catalogue.                           ',
        '                                                ',
        'apply_up_down_stream_correction    False        ',
        '# True/False                                    ',
        'exclude_initial_orbits             1            ',
        '# number                                        ',
        'exclude_final_orbits               0            ',
        '# number                                        ',
        'exclude_initial_orbit_points       0            ',
        '# number                                        ',
        '                                                ',
        'mcmc_iterations                    300000       ',
        '# number                                        ',
        'mcmc_walkers                       50           ',
        '# number                                        ',
        'mcmc_burned_iterations             200000       ',
        '# number                                        ',
        'spectral_mcmc_iterations           60000        ',
        '# number                                        ',
        'spectral_mcmc_walkers              50           ',
        '# number                                        ',
        'spectral_mcmc_burned_iterations    10000        ',
        '# number                                        ',
        '                                                ',
        'first_orbit_ramp                   True         ',
        '# True/False                                    ',
        'second_order_ramp                  False        ',
        '# True/False                                    ',
        'mid_orbit_ramps                    True         ',
        '# True/False                                    ',
        '                                                ',
        'fit_ldc1                           False        ',
        '# True/False                                    ',
        'fit_ldc2                           False        ',
        '# True/False                                    ',
        'fit_ldc3                           False        ',
        '# True/False                                    ',
        'fit_ldc4                           False        ',
        '# True/False                                    ',
        'fit_inclination                    False        ',
        '# True/False                                    ',
        'fit_sma_over_rs                    False        ',
        '# True/False                                    ',
        'fit_mid_time                       True         '
    ]))
    w.close()

    os.system('iraclis -p {0}'.format(parameters_file))
