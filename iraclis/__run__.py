"""Usage: __run__.py -p PARFILE
          __run__.py -d DATADIR
          __run__.py -D DATADIR [PROCEDURE] [PARSTRING]
          __run__.py -T

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .images1_reduction import *
from .images2_calibration import *
from .images3_photometry import *
from .lightcurves1_fitting import *


def process_visit(parameters_file=None, data_directory=None, procedure=None, par_string=None):

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

    print('')

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
                new_frame = [pf.PrimaryHDU(header=frame[0].header, data=frame[0].data),
                             pf.ImageHDU(header=frame2[1 + sample * 5].header, data=frame2[1 + sample * 5].data),
                             pf.ImageHDU(header=frame2[2 + sample * 5].header, data=frame2[2 + sample * 5].data),
                             pf.ImageHDU(header=frame2[3 + sample * 5].header, data=frame2[3 + sample * 5].data),
                             pf.ImageHDU(header=frame2[4 + sample * 5].header, data=frame2[4 + sample * 5].data),
                             pf.ImageHDU(header=frame2[5 + sample * 5].header, data=frame2[5 + sample * 5].data),
                             pf.ImageHDU(header=frame2[6 + sample * 5].header, data=frame2[6 + sample * 5].data),
                             pf.ImageHDU(header=frame2[7 + sample * 5].header, data=frame2[7 + sample * 5].data),
                             pf.ImageHDU(header=frame2[8 + sample * 5].header, data=frame2[8 + sample * 5].data),
                             pf.ImageHDU(header=frame2[9 + sample * 5].header, data=frame2[9 + sample * 5].data),
                             pf.ImageHDU(header=frame2[10 + sample * 5].header, data=frame2[10 + sample * 5].data),
                             pf.ImageHDU(header=frame['SKYAREA'].header, data=frame['SKYAREA'].data,
                                         name='SKYAREA'),
                             pf.ImageHDU(header=frame['SCANMAP'].header, data=frame['SCANMAP'].data,
                                         name='SCANMAP'),
                             pf.ImageHDU(header=frame['WMAP'].header, data=frame['WMAP'].data, name='WMAP'),
                             pf.ImageHDU(header=frame['NWMAP'].header, data=frame['NWMAP'].data, name='NWMAP'),
                             pf.ImageHDU(header=frame['BPCRMAP'].header, data=frame['BPCRMAP'].data,
                                         name='BPCRMAP'),
                             ]

                data_set.spectroscopic_images[i] = pf.HDUList(new_frame)

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

        print('')

    # fitting the white and the spectral light curves

    if variables.fitting_white.value:

        light_curve = plc.open_dict(lc_file)
        light_curve = fitting(light_curve, fitting_spectrum=False)

        print('Saving fitting results in {} ...'.format(str(os.sep).join(fit_file.split(os.sep)[-2:])))
        plc.save_dict(light_curve, fit_file)

        print('Saving fitting plots in {} ...'.format(os.path.split(fitting_figures_directory)[1]))
        plot_fitting(light_curve, fitting_figures_directory)

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

    print ('Processed directory: ', os.path.split(directory)[1])
    print ('                  @: ', os.path.split(directory)[0])
    print('')


def download_visit(vist):
    pass


def white_global_fit_detrended(list_of_files, output_directory,
                               method=False,
                               iterations=200000, walkers=200, burn=150000, precision=3,
                               fit_ld=False, fit_period=False, fit_sma_over_rs=False, fit_inclination=False):

    datasets = list(list_of_files)
    datasets.sort()
    datasets = [pickle.load(open(ff))['lightcurves']['white'] for ff in datasets]

    if fit_period:
        mid_time = datasets[0]['parameters']['t_0']['value']
        period = datasets[0]['parameters']['P']['value']
        xdata = []
        ydata = []
        yerror = []
        for num, dataset in enumerate(datasets):
            print(dataset['parameters']['t_0']['value'])
            xdata.append(round((dataset['parameters']['t_0']['value'] -
                                datasets[0]['parameters']['t_0']['value'])/period))
            ydata.append(dataset['parameters']['t_0']['value'])
            yerror.append(max(dataset['parameters']['t_0']['m_error'], dataset['parameters']['t_0']['p_error']))

        xdata = np.array(xdata)
        ydata = np.array(ydata)
        yerror = np.array(yerror)

        def model(x, a, b):
            return x * a + b

        best_fit, covariance = curve_fit(model, xdata, ydata, sigma=yerror, p0=[period, mid_time], maxfev=20000)

        period, mid_time = best_fit
        print(period, np.sqrt(covariance[0][0]))

    else:
        period = datasets[0]['parameters']['P']['value']

    data = []
    for num, dataset in enumerate(datasets):
        data.append([(dataset['input_time_series']['hjd'] -
                      dataset['parameters']['t_0']['value'] + 1000.0 + num * period),
                     dataset['input_time_series']['raw_lc']/dataset['output_time_series']['systematics'],
                     dataset['input_time_series']['raw_lc_error']/dataset['output_time_series']['systematics']
                     ])

    datasets = datasets[:1]

    if not method:
        method = datasets[0]['limb_darkening']['method']

    limb_darkening_coefficients = 'fit'
    if not fit_ld:
        if method == 'linear':
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value']]
        elif method in ['quad', 'sqrt']:
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value'],
                                           datasets[0]['parameters']['ldc_2']['value']]
        elif method == 'claret':
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value'],
                                           datasets[0]['parameters']['ldc_2']['value'],
                                           datasets[0]['parameters']['ldc_3']['value'],
                                           datasets[0]['parameters']['ldc_4']['value']]

    fit_rp_over_rs = [datasets[0]['parameters']['rp']['value'] / 2, datasets[0]['parameters']['rp']['value'] * 2]
    if fit_sma_over_rs:
        fit_sma_over_rs = [datasets[0]['parameters']['a']['value'] / 2, datasets[0]['parameters']['a']['value'] * 2]
    if fit_inclination:
        fit_inclination = [datasets[0]['parameters']['i']['value'] - 10, 90]

    mcmc = plc.TransitAndPolyFitting(
        data=data,
        method=method,
        limb_darkening_coefficients=limb_darkening_coefficients,
        rp_over_rs=datasets[0]['parameters']['rp']['value'],
        period=period,
        sma_over_rs=datasets[0]['parameters']['a']['value'],
        eccentricity=datasets[0]['parameters']['e']['value'],
        inclination=datasets[0]['parameters']['i']['value'],
        periastron=datasets[0]['parameters']['omega']['value'],
        mid_time=1000.0,
        iterations=iterations,
        walkers=walkers,
        burn=burn,
        precision=precision,
        time_factor=int(round(datasets[0]['exposure']['exp_time']/datasets[0]['exposure']['model_resolution'])),
        exp_time=datasets[0]['exposure']['exp_time'],
        fit_first_order=False,
        fit_second_order=False,
        fit_rp_over_rs=fit_rp_over_rs,
        fit_sma_over_rs=fit_sma_over_rs,
        fit_inclination=fit_inclination,
        )

    mcmc.run_mcmc()

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    mcmc.save_all(os.path.join(output_directory, 'data_base.pickle'))
    mcmc.save_results(os.path.join(output_directory, 'results.txt'))
    mcmc.plot_corner(os.path.join(output_directory, 'correlations.pdf'))
    mcmc.plot_traces(os.path.join(output_directory, 'traces.pdf'))
    mcmc.plot_models(os.path.join(output_directory, 'full_models.pdf'))
    mcmc.plot_detrended_models(os.path.join(output_directory, 'detrended_models.pdf'))


def spectral_global_fit_detrended(list_of_files, lc_id, output_directory,
                                  method=False,
                                  iterations=200000, walkers=200, burn=150000, precision=3,
                                  fit_ld=False):

    datasets = list(list_of_files)
    datasets.sort()
    datasets = [pickle.load(open(ff, 'rb'))['lightcurves'][lc_id] for ff in datasets]

    data = []
    for dataset in datasets:
        data.append([dataset['input_time_series']['hjd'] - dataset['parameters']['t_0']['value'] + 1000.0,
                     dataset['input_time_series']['relative_lc']/dataset['output_time_series']['systematics'],
                     dataset['input_time_series']['relative_lc_error']/dataset['output_time_series']['systematics']
                     ])

    datasets = datasets[:1]

    if not method:
        method = datasets[0]['limb_darkening']['method']

    limb_darkening_coefficients = 'fit'
    if not fit_ld:
        if method == 'linear':
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value']]
        elif method in ['quad', 'sqrt']:
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value'],
                                           datasets[0]['parameters']['ldc_2']['value']]
        elif method == 'claret':
            limb_darkening_coefficients = [datasets[0]['parameters']['ldc_1']['value'],
                                           datasets[0]['parameters']['ldc_2']['value'],
                                           datasets[0]['parameters']['ldc_3']['value'],
                                           datasets[0]['parameters']['ldc_4']['value']]

    mcmc = plc.TransitAndPolyFitting(
        data=data,
        method=method,
        limb_darkening_coefficients=limb_darkening_coefficients,
        rp_over_rs=datasets[0]['parameters']['rp']['value'],
        period=datasets[0]['parameters']['P']['value'],
        sma_over_rs=datasets[0]['parameters']['a']['value'],
        eccentricity=datasets[0]['parameters']['e']['value'],
        inclination=datasets[0]['parameters']['i']['value'],
        periastron=datasets[0]['parameters']['omega']['value'],
        mid_time=1000.0,
        iterations=iterations,
        walkers=walkers,
        burn=burn,
        precision=precision,
        time_factor=int(round(datasets[0]['exposure']['exp_time']/datasets[0]['exposure']['model_resolution'])),
        exp_time=datasets[0]['exposure']['exp_time'],
        fit_first_order=False,
        fit_second_order=False,
        fit_rp_over_rs=[datasets[0]['parameters']['rp']['value'] / 2, datasets[0]['parameters']['rp']['value'] * 2]
        )

    mcmc.run_mcmc()

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    mcmc.save_all(os.path.join(output_directory, 'data_base.pickle'))
    mcmc.save_results(os.path.join(output_directory, 'results.txt'))
    mcmc.plot_corner(os.path.join(output_directory, 'correlations.pdf'))
    mcmc.plot_traces(os.path.join(output_directory, 'traces.pdf'))
    mcmc.plot_models(os.path.join(output_directory, 'full_models.pdf'))
    mcmc.plot_detrended_models(os.path.join(output_directory, 'detrended_models.pdf'))


def console():

    arguments = docopt.docopt(__doc__)

    if arguments['-p'] or arguments['-d'] or arguments['-D']:
        process_visit(arguments['PARFILE'], arguments['DATADIR'], arguments['PROCEDURE'], arguments['PARSTRING'])

    # for developers use only

    elif arguments['-T']:
        from .__test__ import run_test
        run_test()
