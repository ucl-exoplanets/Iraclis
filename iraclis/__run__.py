"""Usage: __run__.py -p PARFILE
          __run__.py -d DATADIR [PROCEDURE] [RESOLUTION] [OUTPUT]

"""

from images1_reduction import *
from images2_calibration import *
from images3_photometry import *
from lightcurves1_fitting import *


def process_visit(parameters_file=None):

    if not parameters_file:
        arguments = docopt.docopt(__doc__)
        # print(arguments)
        # exit()
        if arguments['-p']:
            parameters_file = arguments['PARFILE']
        pipeline_variables.from_parameters_file(parameters_file)
        if arguments['-d']:
            if arguments['DATADIR']:
                pipeline_variables.data_directory.set(arguments['DATADIR'])
            if arguments['PROCEDURE']:
                pipeline_variables.reduction.set(bool(int(arguments['PROCEDURE'][0])))
                pipeline_variables.splitting.set(bool(int(arguments['PROCEDURE'][1])))
                pipeline_variables.extraction.set(bool(int(arguments['PROCEDURE'][2])))
                pipeline_variables.splitting_extraction.set(bool(int(arguments['PROCEDURE'][3])))
                pipeline_variables.fitting_white.set(bool(int(arguments['PROCEDURE'][4])))
                pipeline_variables.fitting_spectrum.set(bool(int(arguments['PROCEDURE'][5])))
            if arguments['RESOLUTION']:
                if arguments['RESOLUTION'] == 'high':
                    pipeline_variables.bins_file.set('default_high')
                elif arguments['RESOLUTION'] == 'low':
                    pipeline_variables.bins_file.set('default_low')
                elif arguments['RESOLUTION'] == 'vlow':
                    pipeline_variables.bins_file.set('default_vlow')
                else:
                    raise PYWFC3InputError('Wrong default resolution given. Options are: high, low, vlow.')
            if arguments['OUTPUT']:
                pipeline_variables.output_directory_copy.set(arguments['OUTPUT'])

    else:
        pipeline_variables.from_parameters_file(parameters_file)

    # load the parameters file

    directory = os.path.abspath(pipeline_variables.data_directory.value)
    print(pipeline_variables.reduction.value,
          pipeline_variables.splitting.value,
          pipeline_variables.extraction.value,
          pipeline_variables.splitting_extraction.value,
          pipeline_variables.fitting_white.value,
          pipeline_variables.fitting_spectrum.value,
          pipeline_variables.bins_file.value,
          pipeline_variables.output_directory_copy.value)

    if not os.path.isdir(directory):
        raise PYWFC3FileError('No such directory: ' + directory)
    else:
        print ''
        print 'Processing directory: ', os.path.split(directory)[1]
        print '                   @: ', os.path.split(directory)[0]
        print ''

    # set paths

    raw_data_directory = os.path.join(directory, pipeline_variables.raw_data_directory.value)
    reduced_data_directory = os.path.join(directory, pipeline_variables.reduced_data_directory.value)
    output_directory = os.path.join(directory, pipeline_variables.output_directory.value)
    output_directory_copy = os.path.join(directory, pipeline_variables.output_directory_copy.value)
    figures_directory = pipeline_variables.figures_directory.value
    extraction_figures_directory = os.path.join(output_directory, figures_directory + '_extraction')
    fitting_figures_directory = os.path.join(output_directory, figures_directory + '_fitting')
    lc_file = os.path.join(output_directory, pipeline_variables.light_curve_file.value + '.pickle')
    fit_file = os.path.join(output_directory, pipeline_variables.fitting_file.value + '.pickle')
    splitted_data_directory = os.path.join(directory, pipeline_variables.splitted_data_directory.value)
    splitted_reduced_data_directory = os.path.join(splitted_data_directory,
                                                   pipeline_variables.reduced_data_directory.value)

    # set process steps and detect inconsistencies

    if pipeline_variables.fitting_spectrum.value and not pipeline_variables.fitting_white.value:
        if not os.path.isfile(fit_file):
            print '\nResetting Fitting-white to True ...'
            pipeline_variables.fitting_white.set(True)

    if (pipeline_variables.fitting_white.value and not pipeline_variables.extraction.value
       and not pipeline_variables.splitting_extraction.value):
        if not os.path.isfile(lc_file):
            print '\nResetting Extraction to True ...'
            pipeline_variables.extraction.set(True)

    if pipeline_variables.extraction.value and not pipeline_variables.reduction.value:
        if not os.path.isdir(reduced_data_directory):
            print '\nResetting Reduction to True ...'
            pipeline_variables.reduction.set(True)

    if pipeline_variables.splitting_extraction.value and not pipeline_variables.splitting.value:
        if not os.path.isdir(splitted_data_directory):
            print '\nResetting Splitting to True ...'
            pipeline_variables.splitting.set(True)

    if pipeline_variables.splitting_extraction.value and pipeline_variables.extraction.value:
            print '\nResetting Splitting extraction to False ...'
            pipeline_variables.splitting_extraction.set(False)

    if pipeline_variables.splitting.value and not pipeline_variables.reduction.value:
        if not os.path.isdir(reduced_data_directory):
            print '\nResetting Reduction to True ...'
            pipeline_variables.reduction.set(True)

    if pipeline_variables.reduction.value:
        if not os.path.isdir(raw_data_directory):
            print '\nMoving raw images in {} ...'.format(pipeline_variables.raw_data_directory.value)
            os.mkdir(raw_data_directory)
            for fits_file in glob.glob(os.path.join(directory, '*.fits')):
                shutil.move(fits_file, os.path.join(raw_data_directory, os.path.split(fits_file)[1]))

    print ''

    # reduction and calibration

    if pipeline_variables.reduction.value:

        # load raw images

        print 'Loading raw images ...'
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

        print 'Saving reduced frames in {} ...'.format(os.path.split(reduced_data_directory)[1])
        data_set.save(reduced_data_directory)
        del data_set
        print ''

    # splitting

    if pipeline_variables.splitting.value:

        print 'Loading raw images ...'
        raw_data_set = DataSet(raw_data_directory)

        # load reduced images

        print 'Loading reduced images ...'
        data_set = DataSet(reduced_data_directory)

        # split images

        if os.path.isdir(splitted_data_directory):
            backup = splitted_data_directory + '_' + time.strftime('%y-%m-%d_%H-%M-%S')
            shutil.copytree(splitted_data_directory, backup)
            shutil.rmtree(splitted_data_directory)

        os.mkdir(splitted_data_directory)

        total_samples = raw_data_set.spectroscopic_images[0][0].header[pipeline_variables.total_samples.keyword] - 1

        for sample in range(total_samples):

            print 'Splitting sample {0}:'.format(sample + 1)

            for i in range(len(raw_data_set.spectroscopic_images)):

                frame = tools.fits_like(data_set.spectroscopic_images[i])
                frame2 = tools.fits_like(raw_data_set.spectroscopic_images[i])
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
            data_set.save(splitted_reduced_data_directory + '_' + str(sample + 1).zfill(2))

        del data_set
        del raw_data_set

        print ''

    # create a backup if output directory exists

    if (pipeline_variables.extraction.value or pipeline_variables.splitting_extraction.value or
            pipeline_variables.fitting_white.value or pipeline_variables.fitting_spectrum.value):

        if os.path.isdir(output_directory):
            backup = output_directory + '_' + time.strftime('%y-%m-%d_%H-%M-%S')
            shutil.copytree(output_directory, backup)
        else:
            os.mkdir(output_directory)

    # extraction of light curves

    if pipeline_variables.extraction.value or pipeline_variables.splitting_extraction.value:

        if pipeline_variables.extraction.value:

            # load reduced images

            print 'Loading reduced images ...'
            data_set = DataSet(reduced_data_directory)

            # use the photometry function to extract the light curves

            data_set, light_curve = photometry(data_set)

        else:

            print 'Loading reduced images ...'
            data_set = DataSet(reduced_data_directory)

            # use the photometry function to extract the light curves

            data_set, light_curve = photometry(data_set)

            star_y_position_array = pipeline_variables.y_star_array.custom_from_dictionary(light_curve)
            y_shift_error_array = pipeline_variables.y_shift_error_array.custom_from_dictionary(light_curve)
            scan_length_array = pipeline_variables.scan_length_array.custom_from_dictionary(light_curve)
            scan_length_error_array = pipeline_variables.scan_length_error_array.custom_from_dictionary(light_curve)

            # load reduced splitted images

            print 'Loading splitted reduced images ...'
            data_set = DataSet(splitted_data_directory)

            # use the split_photometry function to extract the light curves

            data_set, light_curve = split_photometry(data_set)

            star_y_position_array.to_dictionary(light_curve)
            y_shift_error_array.to_dictionary(light_curve)
            scan_length_array.to_dictionary(light_curve)
            scan_length_error_array.to_dictionary(light_curve)

        # save extraction results

        print 'Saving extracted light-curves in {} ...'.format(str(os.sep).join(lc_file.split(os.sep)[-2:]))
        save_lightcurve(light_curve, lc_file)

        # plot extraction diagnostics and results

        print 'Saving extraction plots in {} ...'.format(os.path.split(extraction_figures_directory)[1])
        plot_photometry(data_set, light_curve, extraction_figures_directory)

        del data_set

        print ''

    # fitting the white and the spectral light curves

    if pipeline_variables.fitting_white.value:

        light_curve = open_lightcurve(lc_file)
        light_curve = fitting(light_curve, fitting_spectrum=False)

        print 'Saving fitting results in {} ...'.format(str(os.sep).join(fit_file.split(os.sep)[-2:]))
        save_lightcurve(light_curve, fit_file)

        print 'Saving fitting plots in {} ...'.format(os.path.split(fitting_figures_directory)[1])
        plot_fitting(light_curve, fitting_figures_directory)

    # fitting the spectral light curves

    if pipeline_variables.fitting_spectrum.value:

        light_curve = open_lightcurve(lc_file)
        fitted_white_light_curve = open_lightcurve(fit_file)
        light_curve = fitting(light_curve, fitted_white_light_curve=fitted_white_light_curve,
                              fitting_spectrum=pipeline_variables.fitting_spectrum.value)

        print 'Saving fitting results in {} ...'.format(str(os.sep).join(fit_file.split(os.sep)[-2:]))
        save_lightcurve(light_curve, fit_file)

        print 'Saving fitting plots in {} ...'.format(os.path.split(fitting_figures_directory)[1])
        plot_fitting(light_curve, fitting_figures_directory)

        print ''

    # copy results
    # create a backup if output directory exists

    if pipeline_variables.output_directory_copy.value != 'False':

        if output_directory_copy != output_directory:

            if os.path.isdir(output_directory_copy):
                backup = output_directory_copy + '_' + time.strftime('%y-%m-%d_%H-%M-%S')
                shutil.copytree(output_directory_copy, backup)
                shutil.rmtree(output_directory_copy)

            print 'Copying {} to {} ...'.format(os.path.split(output_directory)[1],
                                                os.path.split(output_directory_copy)[1])
            if os.path.isdir(output_directory):
                shutil.copytree(output_directory, output_directory_copy)

            print ''

    print 'Processed directory: ', os.path.split(directory)[1]
    print '                  @: ', os.path.split(directory)[0]
    print ''


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
            print dataset['parameters']['t_0']['value']
            xdata.append(round((dataset['parameters']['t_0']['value'] - datasets[0]['parameters']['t_0']['value'])/period))
            ydata.append(dataset['parameters']['t_0']['value'])
            yerror.append(max(dataset['parameters']['t_0']['m_error'], dataset['parameters']['t_0']['p_error']))

        xdata = np.array(xdata)
        ydata = np.array(ydata)
        yerror = np.array(yerror)

        def model(x, a, b):
            return x * a + b

        best_fit, covariance = curve_fit(model, xdata, ydata, sigma=yerror, p0=[period, mid_time], maxfev=20000)

        period, mid_time = best_fit
        print period, np.sqrt(covariance[0][0])

    else:
        period = datasets[0]['parameters']['P']['value']

    data = []
    for num, dataset in enumerate(datasets):
        data.append([dataset['input_time_series']['hjd'] - dataset['parameters']['t_0']['value'] + 1000.0 + num * period,
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
    datasets = [pickle.load(open(ff))['lightcurves'][lc_id] for ff in datasets]

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


if __name__ == '__main__':
    process_visit()
