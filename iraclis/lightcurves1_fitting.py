"""
lightcurves1_fitting.py

Includes all the function that performs the fitting. This takes as input a dictionary that contains the extracted
light-curves (produced by the photometry functions in the images3_extraction.py file) and returns a dictionary
that contains the fitting results.

Functions included:
fitting:    ...

"""

from ._3objects import *


def fitting(light_curve, fitted_white_light_curve=None, fitting_spectrum=True,
            apply_up_down_stream_correction=None,
            exclude_initial_orbits=None, exclude_final_orbits=None, exclude_initial_orbit_points=None,
            mcmc_iterations=None, mcmc_walkers=None, mcmc_burned_iterations=None,
            spectral_mcmc_iterations=None, spectral_mcmc_walkers=None, spectral_mcmc_burned_iterations=None,
            first_orbit_ramp=None, second_order_ramp=None, mid_orbit_ramps=None,
            planet=None, method=None, star_teff=None, star_logg=None, star_meta=None,
            white_ldc1=None, white_ldc2=None, white_ldc3=None, white_ldc4=None, bins_file=None,
            fit_ldc1=None, fit_ldc2=None, fit_ldc3=None, fit_ldc4=None,
            rp_over_rs=None, fp_over_fs=None, period=None,
            sma_over_rs=None, fit_sma_over_rs=None, eccentricity=None, inclination=None, fit_inclination=None,
            periastron=None, mid_time=None, fit_mid_time=None):

    # set pipeline variables from input keywords

    apply_up_down_stream_correction = variables.apply_up_down_stream_correction.custom(
        apply_up_down_stream_correction)
    exclude_initial_orbits = variables.exclude_initial_orbits.custom(exclude_initial_orbits)
    exclude_final_orbits = variables.exclude_final_orbits.custom(exclude_final_orbits)
    exclude_initial_orbit_points = variables.exclude_initial_orbit_points.custom(exclude_initial_orbit_points)
    mcmc_iterations = variables.mcmc_iterations.custom(mcmc_iterations)
    mcmc_walkers = variables.mcmc_walkers.custom(mcmc_walkers)
    mcmc_burned_iterations = variables.mcmc_burned_iterations.custom(mcmc_burned_iterations)
    spectral_mcmc_iterations = variables.spectral_mcmc_iterations.custom(spectral_mcmc_iterations)
    spectral_mcmc_walkers = variables.spectral_mcmc_walkers.custom(spectral_mcmc_walkers)
    spectral_mcmc_burned_iterations = variables.spectral_mcmc_burned_iterations.custom(
        spectral_mcmc_burned_iterations)
    first_orbit_ramp = variables.first_orbit_ramp.custom(first_orbit_ramp)
    second_order_ramp = variables.second_order_ramp.custom(second_order_ramp)
    mid_orbit_ramps = variables.mid_orbit_ramps.custom(mid_orbit_ramps)
    planet = variables.planet.custom(planet)
    method = variables.method.custom(method)
    star_teff = variables.star_teff.custom(star_teff)
    star_logg = variables.star_logg.custom(star_logg)
    star_meta = variables.star_meta.custom(star_meta)
    white_ldc1 = variables.white_ldc1.custom(white_ldc1)
    white_ldc2 = variables.white_ldc2.custom(white_ldc2)
    white_ldc3 = variables.white_ldc3.custom(white_ldc3)
    white_ldc4 = variables.white_ldc4.custom(white_ldc4)
    bins_file = variables.bins_file.custom(bins_file)
    fit_ldc1 = variables.fit_ldc1.custom(fit_ldc1)
    fit_ldc2 = variables.fit_ldc2.custom(fit_ldc2)
    fit_ldc3 = variables.fit_ldc3.custom(fit_ldc3)
    fit_ldc4 = variables.fit_ldc4.custom(fit_ldc4)
    rp_over_rs = variables.rp_over_rs.custom(rp_over_rs)
    fp_over_fs = variables.fp_over_fs.custom(fp_over_fs)
    period = variables.period.custom(period)
    sma_over_rs = variables.sma_over_rs.custom(sma_over_rs)
    fit_sma_over_rs = variables.fit_sma_over_rs.custom(fit_sma_over_rs)
    eccentricity = variables.eccentricity.custom(eccentricity)
    inclination = variables.inclination.custom(inclination)
    fit_inclination = variables.fit_inclination.custom(fit_inclination)
    periastron = variables.periastron.custom(periastron)
    mid_time = variables.mid_time.custom(mid_time)
    fit_mid_time = variables.fit_mid_time.custom(fit_mid_time)

    # load pipeline variables to be used

    ra_target = variables.ra_target.custom()
    ra_target.from_dictionary(light_curve)

    dec_target = variables.dec_target.custom()
    dec_target.from_dictionary(light_curve)

    exposure_time = variables.exposure_time.custom()
    exposure_time.from_dictionary(light_curve)

    bins_number = variables.bins_number.custom()
    bins_number.from_dictionary(light_curve)

    subarray_size = variables.sub_array_size.custom()
    subarray_size.from_dictionary(light_curve)

    grism = variables.grism.custom()
    grism.from_dictionary(light_curve)

    heliocentric_julian_date_array = variables.heliocentric_julian_date_array.custom()
    heliocentric_julian_date_array.from_dictionary(light_curve)

    spectrum_direction_array = variables.spectrum_direction_array.custom()
    spectrum_direction_array.from_dictionary(light_curve)

    sky_background_level_array = variables.sky_background_level_array.custom()
    sky_background_level_array.from_dictionary(light_curve)

    star_x_position_array = variables.x_star_array.custom()
    star_x_position_array.from_dictionary(light_curve)

    x_shift_error_array = variables.x_shift_error_array.custom()
    x_shift_error_array.from_dictionary(light_curve)

    star_y_position_array = variables.y_star_array.custom()
    star_y_position_array.from_dictionary(light_curve)

    y_shift_error_array = variables.y_shift_error_array.custom()
    y_shift_error_array.from_dictionary(light_curve)

    scan_length_array = variables.scan_length_array.custom()
    scan_length_array.from_dictionary(light_curve)

    scan_length_error_array = variables.scan_length_error_array.custom()
    scan_length_error_array.from_dictionary(light_curve)

    white_lower_wavelength = variables.white_upper_wavelength
    white_upper_wavelength = variables.white_upper_wavelength

    flux_array = variables.flux_array.custom()
    error_array = variables.error_array.custom()
    ph_error_array = variables.ph_error_array.custom()
    lower_wavelength = variables.lower_wavelength.custom()
    upper_wavelength = variables.upper_wavelength.custom()
    ldc1 = variables.ldc1.custom()
    ldc2 = variables.ldc2.custom()
    ldc3 = variables.ldc3.custom()
    ldc4 = variables.ldc4.custom()

    # set bins

    white_dictionary, bins_dictionaries = \
        variables.set_binning(light_curve, white_lower_wavelength.value, white_upper_wavelength.value,
                              white_ldc1.value, white_ldc2.value, white_ldc3.value, white_ldc4.value,
                              bins_file.value)

    white_dictionary_x, bins_dictionaries_x = \
        variables.set_binning(light_curve, white_lower_wavelength.value, white_upper_wavelength.value,
                              white_ldc1.value, white_ldc2.value, white_ldc3.value, white_ldc4.value,
                              bins_file.value)

    if bins_number.value != len(bins_dictionaries):
        raise IraclisFileError('Bin file used does not match the extracted bins.')

    for i, j in enumerate([white_dictionary] + bins_dictionaries):
        j.from_dictionary(light_curve)
        ldc1.from_dictionary(([white_dictionary_x] + bins_dictionaries_x)[i])
        ldc1.to_dictionary(j)
        ldc2.from_dictionary(([white_dictionary_x] + bins_dictionaries_x)[i])
        ldc2.to_dictionary(j)
        ldc3.from_dictionary(([white_dictionary_x] + bins_dictionaries_x)[i])
        ldc3.to_dictionary(j)
        ldc4.from_dictionary(([white_dictionary_x] + bins_dictionaries_x)[i])
        ldc4.to_dictionary(j)

    # up-stream / down-stream correction

    if apply_up_down_stream_correction.value:

        test1 = star_y_position_array.value[0] - 507
        test2 = test1 + spectrum_direction_array.value[0] * scan_length_array.value[0]

        if test1 * test2 < 0:
            apply_up_down_stream_correction.set(True)
            print('Correcting for up-stream/down-stream effect ... ')
        else:
            apply_up_down_stream_correction.set(False)
            print('No correction for up-stream/down-stream effect is needed ... ')

    if apply_up_down_stream_correction.value:

        # for scan_direction in [1.0, -1.0]:
        #
        #     fr = np.where(spectrum_direction_array.value == scan_direction)[0]
        #
        #     if len(fr) > 0:
        #
        #         zerofr = star_y_position_array.value[fr]
        #         sigmfr = scan_length_array.value[fr]
        #         begnfr = zerofr
        #         fitfr = np.poly1d(np.polyfit(begnfr, sigmfr, 1))
        #
        #         # if scan_direction > 0:
        #         #     ur = (512.0 / 261.0) * 512.0 / 2.93
        #         #     e = 8.774724
        #         #     Y = 522.0
        #         #
        #         #     l1 = fitfr(begnfr)
        #         #     y1 = begnfr - (1014.0 - 522.0) / 2
        #         #
        #         #     correction = np.nanmean(((l1[0] / (- (l1[0] - Y + 2 * y1[0]) / ur + e)) / (
        #               l1 / (- (l1 - Y + 2 * y1) / ur + e)) - 1.0) / (y1 - y1[0]))
        #         #     print correction
        #         #     correction = 0.000095
        #
        #         for i in [white_dictionary] + bins_dictionaries:
        #
        #             detector_length = subarray_size.value
        #
        #             ur = (512.0 / (detector_length / 2.0)) * 512.0 / 2.93
        #             e = 8.774724
        #             l1 = fitfr(begnfr)
        #             if scan_direction > 0:
        #                 y1 = begnfr - (1014.0 - detector_length) / 2
        #             else:
        #                 y1 = detector_length - (begnfr - (1014.0 - detector_length) / 2)
        #
        #             flux_array.from_dictionary(i)
        #             flux = flux_array.value
        #             # flux[fr] = flux[fr] / ((fitfr(begnfr) / fitfr(begnfr[0])) * (
        #             1 + correction * (begnfr - begnfr[0])))
        #             flux[fr] = flux[fr] / (- (l1 - detector_length + 2 * y1) / ur + e)
        #             # plt.plot(flux, 'o')
        #             # plt.show()
        #             flux_array.set(flux)
        #             flux_array.to_dictionary(i)
        #
        #             error_array.from_dictionary(i)
        #             error = error_array.value
        #             # error[fr] = error[fr] / ((fitfr(begnfr) / fitfr(begnfr[0])) * (
        #             1 + correction * (begnfr - begnfr[0])))
        #             error[fr] = error[fr] / (- (l1 - detector_length + 2 * y1) / ur + e)
        #             error_array.set(error)
        #             error_array.to_dictionary(i)

        for scan_direction in [1.0, -1.0]:

            fr = np.where(spectrum_direction_array.value == scan_direction)[0]

            if len(fr) > 0:

                zerofr = star_y_position_array.value[fr]
                sigmfr = scan_length_array.value[fr]
                begnfr = zerofr
                fitfr = np.poly1d(np.polyfit(begnfr, sigmfr, 1))

                for i in [white_dictionary] + bins_dictionaries:

                    flux_array.from_dictionary(i)
                    flux = flux_array.value
                    flux[fr] = flux[fr] * fitfr(begnfr[0]) / fitfr(begnfr)
                    flux_array.set(flux)
                    flux_array.to_dictionary(i)

                    error_array.from_dictionary(i)
                    error = error_array.value
                    error[fr] = error[fr] * fitfr(begnfr[0]) / fitfr(begnfr)
                    error_array.set(error)
                    error_array.to_dictionary(i)

                    ph_error_array.from_dictionary(i)
                    ph_error = ph_error_array.value
                    ph_error[fr] = ph_error[fr] * fitfr(begnfr[0]) / fitfr(begnfr)
                    ph_error_array.set(ph_error)
                    ph_error_array.to_dictionary(i)

    # exclude orbits / points

    indices_to_remain = np.arange(len(heliocentric_julian_date_array.value))

    if exclude_initial_orbits.value > 0:
        print('Excluding {0} orbit{1} from the beginning of the visit ...'.format(
            exclude_initial_orbits.value, ['s', ''][1 // exclude_initial_orbits.value]))
        htime = heliocentric_julian_date_array.value
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]

        indices_to_remain = indices_to_remain[orbits[exclude_initial_orbits.value]:]

    if exclude_final_orbits.value > 0:
        print('Excluding {0} orbit{1} from the end of the visit ...'.format(
            exclude_final_orbits.value, ['s', ''][1 // exclude_final_orbits.value]))
        htime = heliocentric_julian_date_array.value[indices_to_remain]
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]

        indices_to_remain = indices_to_remain[:orbits[-exclude_final_orbits.value]]

    if exclude_initial_orbit_points.value > 0:
        print('Excluding {0} point{1} from the beginning of each orbit ...'.format(
            exclude_initial_orbit_points.value, ['s', ''][1 // exclude_initial_orbit_points.value]))
        htime = heliocentric_julian_date_array.value[indices_to_remain]
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]

        indices_to_remain = np.delete(indices_to_remain,
                                      np.concatenate([orbits + i for i in range(exclude_initial_orbit_points.value)]))

    for i in [heliocentric_julian_date_array, spectrum_direction_array, sky_background_level_array,
              star_x_position_array, x_shift_error_array, star_y_position_array, y_shift_error_array, scan_length_array,
              scan_length_error_array]:

        i.set(i.value[indices_to_remain])
        i.to_dictionary(light_curve)

    for i in [white_dictionary] + bins_dictionaries:

        flux_array.from_dictionary(i)
        flux_array.set(flux_array.value[indices_to_remain])
        flux_array.to_dictionary(i)

        error_array.from_dictionary(i)
        error_array.set(error_array.value[indices_to_remain])
        error_array.to_dictionary(i)

        ph_error_array.from_dictionary(i)
        ph_error_array.set(ph_error_array.value[indices_to_remain])
        ph_error_array.to_dictionary(i)

    # define hst orbital phases

    if mid_orbit_ramps.value:
        htime = heliocentric_julian_date_array.value
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
        dumps = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
        dphase = np.zeros(len(htime))
        for i in range(1, len(dumps)):
            if dumps[i] not in orbits:
                if i != len(dumps) - 1:
                    for j in range(dumps[i], dumps[i + 1]):
                        dphase[j] = 1
                else:
                    for j in range(dumps[i], len(dphase)):
                        dphase[j] = 1
    else:
        htime = heliocentric_julian_date_array.value
        dphase = np.zeros(len(htime))

    if first_orbit_ramp.value:
        htime = heliocentric_julian_date_array.value
        if mid_orbit_ramps.value:
            orbits = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
        else:
            orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
        orbits = htime[orbits]
        fphase = np.where(htime < orbits[1], 1, 0)

    else:
        htime = heliocentric_julian_date_array.value
        fphase = np.zeros(len(htime))

    htime = heliocentric_julian_date_array.value
    if mid_orbit_ramps.value:
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 5.0 / 60.0 / 24.0)[0]
    else:
        orbits = np.where(abs(htime - np.roll(htime, 1)) > 30.0 / 60.0 / 24.0)[0]
    t0s = htime[orbits]
    ophase = []
    for pp in t0s:
        ppp = htime - pp
        ppp = np.where(ppp < 0, 1000, ppp)
        ophase.append(ppp)

    ophase = np.min(ophase, 0)

    # outliers filter

    if fitted_white_light_curve is None:

        lightcurves = []
        for i in bins_dictionaries:
            flux_array.from_dictionary(i)
            lightcurves.append(flux_array.value)

        ica = FastICA(n_components=len(lightcurves), max_iter=10000)
        components = ica.fit_transform(np.array(lightcurves).T).T

        indices_to_remain = []
        for i in components:
            indices_to_remain.append(np.array(np.abs(i - np.median(i)) < 20 * np.median(np.abs(i - np.median(i)))))
        indices_to_remain = np.where(np.product(indices_to_remain, 0))[0]
        indices_to_remain = np.sort(np.unique(np.array(indices_to_remain)))

        for i in [heliocentric_julian_date_array, spectrum_direction_array, sky_background_level_array,
                  star_x_position_array, x_shift_error_array, star_y_position_array, y_shift_error_array,
                  scan_length_array, scan_length_error_array]:

            i.set(i.value[indices_to_remain])
            i.to_dictionary(light_curve)

        for i in [white_dictionary] + bins_dictionaries:

            flux_array.from_dictionary(i)
            flux_array.set(flux_array.value[indices_to_remain])
            flux_array.to_dictionary(i)

            error_array.from_dictionary(i)
            error_array.set(error_array.value[indices_to_remain])
            error_array.to_dictionary(i)

            ph_error_array.from_dictionary(i)
            ph_error_array.set(ph_error_array.value[indices_to_remain])
            ph_error_array.to_dictionary(i)

        ophase = ophase[indices_to_remain]
        dphase = dphase[indices_to_remain]
        fphase = fphase[indices_to_remain]

    else:

        white_fit_hjd = fitted_white_light_curve['lightcurves']['white']['input_time_series']['hjd']
        htime = heliocentric_julian_date_array.value
        indices_to_remain = np.where([(ff in white_fit_hjd) for ff in htime])[0]

        for i in [heliocentric_julian_date_array, spectrum_direction_array, sky_background_level_array,
                  star_x_position_array, x_shift_error_array, star_y_position_array, y_shift_error_array,
                  scan_length_array,
                  scan_length_error_array]:
            i.set(i.value[indices_to_remain])
            i.to_dictionary(light_curve)

        for i in [white_dictionary] + bins_dictionaries:
            flux_array.from_dictionary(i)
            flux_array.set(flux_array.value[indices_to_remain])
            flux_array.to_dictionary(i)

            error_array.from_dictionary(i)
            error_array.set(error_array.value[indices_to_remain])
            error_array.to_dictionary(i)

            ph_error_array.from_dictionary(i)
            ph_error_array.set(ph_error_array.value[indices_to_remain])
            ph_error_array.to_dictionary(i)

        ophase = ophase[indices_to_remain]
        dphase = dphase[indices_to_remain]
        fphase = fphase[indices_to_remain]

    # match forward and reverse scans

    fr = np.where(spectrum_direction_array.value > 0)[0]
    if len(fr) != len(spectrum_direction_array.value):

        print('Matching forward and revers scans ...')

        fr_out = np.where(spectrum_direction_array.value > 0)[0]
        rv_out = np.where(spectrum_direction_array.value < 0)[0]

        flux_array.from_dictionary(white_dictionary)
        flux = flux_array.value
        shift = np.mean(flux[fr_out]) / np.mean(flux[rv_out])

        for i in [white_dictionary] + bins_dictionaries:

            flux_array.from_dictionary(i)
            flux = flux_array.value
            flux[fr] = flux[fr] / shift
            flux_array.set(flux)
            flux_array.to_dictionary(i)

            error_array.from_dictionary(i)
            error = error_array.value
            error[fr] = error[fr] / shift
            error_array.set(error)
            error_array.to_dictionary(i)

            ph_error_array.from_dictionary(i)
            ph_error = ph_error_array.value
            ph_error[fr] = ph_error[fr] / shift
            ph_error_array.set(ph_error)
            ph_error_array.to_dictionary(i)

    # set target

    if 'auto' in [rp_over_rs.value, fp_over_fs.value, period.value, sma_over_rs.value, eccentricity.value,
                  inclination.value, periastron.value, mid_time.value]:

        catalogue = plc.oec_catalogue()

        if planet.value == 'auto':

            planets = []

            for i in catalogue.planets:
                if not np.isnan(i.system.dec):
                    planets.append([np.sqrt((i.system.dec.deg - dec_target.value) ** 2
                                            + (i.system.ra.deg - ra_target.value) ** 2), i.name, i])
            planets.sort()

            planet.set(planets[0][2].name)

        else:
            planet.set(planet.value.replace('_', ' '))

        oec_parameters = plc.find_oec_parameters(planet.value, catalogue)

        print('Some parameters are set to auto, resetting to the OEC values for {0} ...'.format(planet.value))

        if star_logg.value == 'auto':
            star_logg.set(oec_parameters[1])

        if star_teff.value == 'auto':
            star_teff.set(oec_parameters[2])

        if star_meta.value == 'auto':
            star_meta.set(oec_parameters[3])

        if rp_over_rs.value == 'auto':
            rp_over_rs.set(oec_parameters[4])

        if fp_over_fs.value == 'auto':
            fp_over_fs.set(oec_parameters[5])

        if period.value == 'auto':
            period.set(oec_parameters[6])

        if sma_over_rs.value == 'auto':
            sma_over_rs.set(oec_parameters[7])

        if eccentricity.value == 'auto':
            eccentricity.set(oec_parameters[8])

        if inclination.value == 'auto':
            inclination.set(oec_parameters[9])

        if periastron.value == 'auto':
            periastron.set(oec_parameters[10])

        if mid_time.value == 'auto':
            mid_time.set(oec_parameters[11])

    for j in range(len(bins_dictionaries)):

        ldc1.from_dictionary(bins_dictionaries[j])

        if isinstance(ldc1.value, str):

            method.set('claret')
            if star_teff.value < 4000.0:
                stellar_model = 'PHOENIX'
            else:
                stellar_model = 'ATLAS'

            if ldc1.value == 'default_high' and grism.value == 'G141':
                band = 'WFC3G141H'
            elif ldc1.value == 'default_low' and grism.value == 'G141':
                band = 'WFC3G141L'
            elif ldc1.value == 'default_vlow' and grism.value == 'G141':
                band = 'WFC3G141VL'
            elif ldc1.value == 'default_high' and grism.value == 'G102':
                band = 'WFC3G102H'
            elif ldc1.value == 'default_low' and grism.value == 'G102':
                band = 'WFC3G102L'
            else:
                band = 'WFC3G102VL'

            a1, a2, a3, a4 = plc.clablimb('claret', star_logg.value, star_teff.value, 0.0, '{0}{1}'.format(band, j),
                                          stellar_model=stellar_model)

            ldc1.to_dictionary(bins_dictionaries[j], value=a1)
            ldc2.to_dictionary(bins_dictionaries[j], value=a2)
            ldc3.to_dictionary(bins_dictionaries[j], value=a3)
            ldc4.to_dictionary(bins_dictionaries[j], value=a4)

    ldc1.from_dictionary(white_dictionary)

    if ldc1.value == 'default':
        method.set('claret')
        if star_teff.value < 4000.0:
            stellar_model = 'PHOENIX'
        else:
            stellar_model = 'ATLAS'

        if grism.value == 'G141':
            band = 'WFC3G141W'
        else:
            band = 'WFC3G102W'

        a1, a2, a3, a4 = plc.clablimb('claret', star_logg.value, star_teff.value, 0.0, band,
                                      stellar_model=stellar_model)

        ldc1.to_dictionary(white_dictionary, value=a1)
        ldc2.to_dictionary(white_dictionary, value=a2)
        ldc3.to_dictionary(white_dictionary, value=a3)
        ldc4.to_dictionary(white_dictionary, value=a4)

    # set observation

    eclipse_time = mid_time.value + 0.5 * period.value
    meanhjd_time = np.mean(heliocentric_julian_date_array.value)
    transit_test = abs(mid_time.value - meanhjd_time) / period.value
    transit_test = abs(transit_test - np.round(transit_test))
    eclipse_test = abs(eclipse_time - meanhjd_time) / period.value
    eclipse_test = abs(eclipse_test - np.round(eclipse_test))

    if transit_test < eclipse_test:
        observation_type = 'transit'
        print('This is a transit observation.')
    else:
        observation_type = 'eclipse'
        print('This is an eclipse observation.')

    # import series

    data_time = heliocentric_julian_date_array.value
    flux_array.from_dictionary(white_dictionary)
    data_white = flux_array.value
    error_array.from_dictionary(white_dictionary)
    data_white_error = error_array.value
    ph_error_array.from_dictionary(white_dictionary)
    data_white_ph_error = ph_error_array.value
    data_ophase = np.array(ophase)
    data_dphase = np.array(dphase)
    data_fphase = np.array(fphase)
    data_scan = spectrum_direction_array.value
    data_sky = sky_background_level_array.value
    data_y_shift = star_y_position_array.value
    data_y_shift_error = y_shift_error_array.value
    data_x_shift = star_x_position_array.value
    data_x_shift_error = x_shift_error_array.value

    names = []
    print_names = []
    limits1 = []
    limits2 = []
    initial = []

    # forward scans normalisation factors
    names.append('n_w_for')
    print_names.append('n_\mathrm{w}^\mathrm{for}')
    initial.append(np.median(data_white))
    if (data_scan < 0).all():
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(np.median(data_white) / 10.0)
        limits2.append(np.median(data_white) * 10.0)

    # reverse scans normalisation factors
    names.append('n_w_rev')
    print_names.append('n_\mathrm{w}^\mathrm{rev}')
    initial.append(np.median(data_white))
    if (data_scan > 0).all():
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(np.median(data_white) / 10.0)
        limits2.append(np.median(data_white) * 10.0)

    # long term ramp - 1st order
    names.append('r_a1')
    print_names.append('r_{a1}')
    initial.append(0.001)
    limits1.append(-10.0)
    limits2.append(10.0)

    # long term ramp - 2nd order
    names.append('r_a2')
    print_names.append('r_{a2}')
    initial.append(0.0)
    if second_order_ramp.value:
        limits1.append(-10.0)
        limits2.append(10.0)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # sort term ramp - amplitude
    names.append('r_b1')
    print_names.append('r_{b1}')
    initial.append(0.001)
    limits1.append(-10.0)
    limits2.append(10.0)

    # sort term mid-orbit ramp - amplitude
    names.append('mor_b1')
    print_names.append('mor_{b1}')
    initial.append(0.001)
    if np.sum(data_dphase ** 2) == 0:
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(-10.0)
        limits2.append(10.0)

    # sort term first-orbit ramp - amplitude
    names.append('for_b1')
    print_names.append('for_{b1}')
    initial.append(0.001)
    if np.sum(data_fphase ** 2) == 0:
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(-10.0)
        limits2.append(10.0)

    # sort term ramp - decay
    names.append('r_b2')
    print_names.append('r_{b2}')
    initial.append(500.0)
    limits1.append(1.0)
    limits2.append(3000.0)

    # sort term mid-orbit ramp - decay
    names.append('mor_b2')
    print_names.append('mor_{b2}')
    initial.append(500.0)
    if np.sum(data_dphase ** 2) == 0:
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(1.0)
        limits2.append(3000.0)

    # sort term first-orbit ramp - decay
    names.append('for_b2')
    print_names.append('for_{b2}')
    initial.append(500.0)
    if np.sum(data_fphase ** 2) == 0:
        limits1.append(np.nan)
        limits2.append(np.nan)
    else:
        limits1.append(1.0)
        limits2.append(3000.0)

    # ld1
    ldc1.from_dictionary(white_dictionary)
    names.append('ldc_1')
    print_names.append('ldc_1')
    initial.append(ldc1.value)
    if fit_ldc1.value:
        limits1.append(0.000001)
        limits2.append(0.999999)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # ld2
    ldc2.from_dictionary(white_dictionary)
    names.append('ldc_2')
    print_names.append('ldc_2')
    initial.append(ldc2.value)
    if fit_ldc2.value:
        limits1.append(0.000001)
        limits2.append(0.999999)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # ld3
    ldc3.from_dictionary(white_dictionary)
    names.append('ldc_3')
    print_names.append('ldc_3')
    initial.append(ldc3.value)
    if fit_ldc3.value:
        limits1.append(0.000001)
        limits2.append(0.999999)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # ld4
    ldc4.from_dictionary(white_dictionary)
    names.append('ldc_4')
    print_names.append('ldc_4')
    initial.append(ldc4.value)
    if fit_ldc4.value:
        limits1.append(0.000001)
        limits2.append(0.999999)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # rp/rs
    names.append('rp')
    print_names.append('R_\mathrm{p}/R_*')
    initial.append(rp_over_rs.value)
    if observation_type == 'transit':
        limits1.append(rp_over_rs.value / 2)
        limits2.append(rp_over_rs.value * 2)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # fp/rs
    names.append('fp')
    print_names.append('F_\mathrm{p}/F_*')
    initial.append(fp_over_fs.value)
    if observation_type == 'eclipse':
        limits1.append(0.0000001)
        limits2.append(0.1)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # period
    names.append('P')
    print_names.append('P')
    initial.append(period.value)
    limits1.append(np.nan)
    limits2.append(np.nan)

    # semi-major axis
    names.append('a')
    print_names.append('a/R_*')
    initial.append(sma_over_rs.value)
    if fit_sma_over_rs.value:
        limits1.append(max(1.0, sma_over_rs.value / 2))
        limits2.append(sma_over_rs.value * 2)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # eccentricity
    names.append('e')
    print_names.append('e')
    initial.append(eccentricity.value)
    limits1.append(np.nan)
    limits2.append(np.nan)

    # inclination
    names.append('i')
    print_names.append('i')
    initial.append(inclination.value)
    if fit_inclination.value:
        limits1.append(60.0)
        limits2.append(90.0)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # periastron
    names.append('omega')
    print_names.append('\omega')
    initial.append(periastron.value)
    limits1.append(np.nan)
    limits2.append(np.nan)

    # mid-transit time
    names.append('t_0')
    print_names.append('T_0')
    mid_time.set(mid_time.value + period.value * round((data_time[-1] - mid_time.value) / period.value))
    if observation_type == 'eclipse':
        mid_time.set(mid_time.value - period.value / 2)
    initial.append(mid_time.value)
    if fit_mid_time.value:
        limits1.append(mid_time.value - 0.2)
        limits2.append(mid_time.value + 0.2)
    else:
        limits1.append(np.nan)
        limits2.append(np.nan)

    # set model
    sub_exposures = int(exposure_time.value / 10.0 + 1)
    sub_exposure_time = float(exposure_time.value / sub_exposures)
    sub_exposure_time /= (24.0 * 60 * 60)

    if exposure_time.value > 10:
        exposure_time.set(exposure_time.value / (24.0 * 60 * 60))
        data_time_hr = np.array([])
        for i in range(sub_exposures):
            data_time_hr = np.append(data_time_hr,
                                     data_time - exposure_time.value / 2.0 + (i + 0.5) * sub_exposure_time)
    else:
        data_time_hr = data_time

    def model(model_norm_f, model_norm_r, model_r_a1, model_r_a2,
              model_r_b1, model_mor_b1, model_for_b1, model_r_b2, model_mor_b2, model_for_b2,
              model_ldc1, model_ldc2, model_ldc3, model_ldc4,
              model_rp_over_rs, model_fp_over_fs, model_period, model_sma_over_rs, model_eccentricity,
              model_inclination, model_periastron, model_mid_time, independent=False):

        normalisation = np.where(data_scan > 0, model_norm_f, model_norm_r)
        detrend1 = (1.0 - model_r_a1 * (data_time - model_mid_time)
                    + model_r_a2 * ((data_time - model_mid_time) ** 2))
        ramp_ampl = np.where(data_dphase == 0, model_r_b1, model_mor_b1)
        ramp_ampl = np.where(data_fphase == 0, ramp_ampl, model_for_b1)
        ramp_decay = np.where(data_dphase == 0, model_r_b2, model_mor_b2)
        ramp_decay = np.where(data_fphase == 0, ramp_decay, model_for_b2)
        detrend2 = 1.0 - ramp_ampl * np.exp(- ramp_decay * data_ophase)

        if observation_type == 'transit':
            signal = plc.transit(method.value, [model_ldc1, model_ldc2, model_ldc3, model_ldc4],
                                 model_rp_over_rs, model_period,
                                 model_sma_over_rs, model_eccentricity, model_inclination,
                                 model_periastron, model_mid_time, data_time_hr)
        else:
            signal = plc.eclipse(model_fp_over_fs, model_rp_over_rs, model_period,
                                 model_sma_over_rs, model_eccentricity, model_inclination,
                                 model_periastron, model_mid_time + model_period / 2, data_time_hr)

        signal = np.mean(np.reshape(signal, (sub_exposures, len(data_time))), 0)

        if not independent:
            return normalisation * detrend1 * detrend2 * signal
        else:
            return normalisation * detrend1 * detrend2, signal

    fitting_results = {'lightcurves': {}, 'spectrum': {}}

    if fitted_white_light_curve is None:

        mcmc = plc.EmceeFitting([data_white, data_white_error], model, np.array(initial), np.array(limits1),
                                np.array(limits2), mcmc_walkers.value, mcmc_iterations.value,
                                mcmc_burned_iterations.value, names, print_names, counter=True)

        mcmc.run_mcmc()

        data_white_error *= np.std(data_white - model(*mcmc.results['parameters_final'])) / np.median(data_white_error)

        mcmc = plc.EmceeFitting([data_white, data_white_error], model, np.array(initial), np.array(limits1),
                                np.array(limits2), mcmc_walkers.value, mcmc_iterations.value,
                                mcmc_burned_iterations.value, names, print_names, counter=True)

        mcmc.run_mcmc()

        white_fit = mcmc.results

        lower_wavelength.from_dictionary(white_dictionary)
        upper_wavelength.from_dictionary(white_dictionary)

        white_fit['wavelength'] = {'lambda1': lower_wavelength.value,
                                   'lambda2': upper_wavelength.value,
                                   'lambda_mean': 0.5 * (lower_wavelength.value + upper_wavelength.value),
                                   'lambda_width': upper_wavelength.value - lower_wavelength.value
                                   }

        white_fit['limb_darkening'] = {'method': method.value}

        white_fit['exposure'] = {'exp_time': exposure_time.value * (24.0 * 60 * 60),
                                 'model_resolution': sub_exposure_time * (24.0 * 60 * 60)
                                 }

        white_fit['input_time_series'] = {'hjd': data_time,
                                          'scan': data_scan,
                                          'sky': data_sky,
                                          'x_shift': data_x_shift,
                                          'x_shift_error': data_x_shift_error,
                                          'y_shift': data_y_shift,
                                          'y_shift_error': data_y_shift_error,
                                          'raw_lc': data_white,
                                          'raw_lc_error': data_white_error,
                                          'raw_lc_ph_error': data_white_ph_error
                                          }

        del white_fit['input_series']

        model_systematics, model_curve = model(*mcmc.results['parameters_final'], independent=True)
        white_fit['output_time_series'] = {'phase': ((data_time - white_fit['parameters']['t_0']['value']) /
                                                     white_fit['parameters']['P']['value']),
                                           'full_model': model_systematics * model_curve,
                                           'full_residuals': data_white - model_systematics * model_curve,
                                           'systematics': model_systematics,
                                           'detrended_lc': data_white / model_systematics,
                                           'transit': model_curve,
                                           'residuals': data_white / model_systematics - model_curve
                                           }

        del white_fit['output_series']

        chis = np.round(
            np.sum((white_fit['output_time_series']['full_residuals'] ** 2) / (
                    white_fit['input_time_series']['raw_lc_error'] ** 2)), 1)
        rchis = np.round(
            np.sum((white_fit['output_time_series']['full_residuals'] ** 2) / (
                    white_fit['input_time_series']['raw_lc_error'] ** 2)) / (
                len(data_time) - len(white_fit['statistics']['corr_matrix'][0])), 2)
        std = int(np.round(np.std(white_fit['output_time_series']['residuals']) * 1000000, 0))
        rstd = np.round(np.std(white_fit['output_time_series']['full_residuals']) / (
            np.mean([np.std(np.random.normal(0, white_fit['input_time_series']['raw_lc_ph_error']))
                     for ff in range(1000)])), 2)
        autocor = np.round(np.max(np.abs(white_fit['statistics']['res_autocorr'][1:3])), 2)
        white_fit['res_statistics'] = {'chisqr': chis,
                                       'rchisqr': rchis,
                                       'stdppm': std,
                                       'rstd': rstd,
                                       'autocor': autocor
                                       }

    else:

        print('Loading previous white fitting results')
        white_fit = fitted_white_light_curve['lightcurves']['white']

    fitting_results['lightcurves']['white'] = white_fit

    # fit relative spectrum

    if fitting_spectrum:

        counter = PipelineCounter('Spectrum', len(bins_dictionaries))

        for j in range(len(bins_dictionaries)):

            ldc1.from_dictionary(bins_dictionaries[j])
            ldc2.from_dictionary(bins_dictionaries[j])
            ldc3.from_dictionary(bins_dictionaries[j])
            ldc4.from_dictionary(bins_dictionaries[j])
            model_white = fitting_results['lightcurves']['white']['output_time_series']['transit']
            flux_array.from_dictionary(white_dictionary)
            data_white = flux_array.value
            model_white = data_white / model_white
            model_white = model_white / np.mean(model_white)
            # error_array.from_dictionary(white_dictionary)
            # data_white_error = error_array.value
            flux_array.from_dictionary(bins_dictionaries[j])
            data_bin = flux_array.value
            error_array.from_dictionary(bins_dictionaries[j])
            data_bin_error = error_array.value
            ph_error_array.from_dictionary(bins_dictionaries[j])
            data_bin_ph_error = ph_error_array.value
            # data_relative = data_bin / data_white
            # data_relative_error = np.sqrt((((1.0 / data_white) ** 2) * (data_bin_error ** 2) +
            #                               ((data_bin / data_white / data_white) ** 2) * (data_white_error ** 2)))

            names = []
            print_names = []
            limits1 = []
            limits2 = []
            initial = []

            # forward scans normalisation factors
            names.append('n_l_for')
            print_names.append('n_\lambda ^\mathrm{for}')
            initial.append(np.median(data_bin))
            if (data_scan < 0).all():
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(np.median(data_bin) / 10.0)
                limits2.append(np.median(data_bin) * 10.0)

            # reverse scans normalisation factors
            names.append('n_l_rev')
            print_names.append('n_\lambda ^\mathrm{rev}')
            initial.append(np.median(data_bin))
            if (data_scan > 0).all():
                limits1.append(np.nan)
                limits2.append(np.nan)
            else:
                limits1.append(np.median(data_bin) / 10.0)
                limits2.append(np.median(data_bin) * 10.0)

            # long term ramp - 1st order
            names.append('r_a1')
            print_names.append('r_{a1}')
            initial.append(0.001)
            limits1.append(-10.0)
            limits2.append(10.0)

            # ld1
            names.append('ldc_1')
            print_names.append('ldc_1')
            initial.append(ldc1.value)
            if fit_ldc1.value:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # ld2
            names.append('ldc_2')
            print_names.append('ldc_2')
            initial.append(ldc2.value)
            if fit_ldc2.value:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # ld3
            names.append('ldc_3')
            print_names.append('ldc_3')
            initial.append(ldc3.value)
            if fit_ldc3.value:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # ld1
            names.append('ldc_4')
            print_names.append('ldc_4')
            initial.append(ldc4.value)
            if fit_ldc4.value:
                limits1.append(0.000001)
                limits2.append(0.999999)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # rp/rs
            names.append('rp')
            print_names.append('R_\mathrm{p}/R_*')
            initial.append(rp_over_rs.value)
            if observation_type == 'transit':
                limits1.append(rp_over_rs.value / 2)
                limits2.append(rp_over_rs.value * 2)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # fp/rs
            names.append('fp')
            print_names.append('F_\mathrm{p}/F_*')
            initial.append(fp_over_fs.value)
            if observation_type == 'eclipse':
                limits1.append(0)
                limits2.append(10 * fp_over_fs.value)
            else:
                limits1.append(np.nan)
                limits2.append(np.nan)

            # period
            names.append('P')
            print_names.append('P')
            initial.append(fitting_results['lightcurves']['white']['parameters']['P']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            # semi-major axis
            names.append('a')
            print_names.append('a/R_*')
            initial.append(fitting_results['lightcurves']['white']['parameters']['a']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            # eccentricity
            names.append('e')
            print_names.append('e')
            initial.append(fitting_results['lightcurves']['white']['parameters']['e']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            # inclination
            names.append('i')
            print_names.append('i')
            initial.append(fitting_results['lightcurves']['white']['parameters']['i']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            # periastron
            names.append('omega')
            print_names.append('\omega')
            initial.append(fitting_results['lightcurves']['white']['parameters']['omega']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            # mid-transit time
            names.append('t_0')
            print_names.append('T_0')
            initial.append(fitting_results['lightcurves']['white']['parameters']['t_0']['value'])
            limits1.append(np.nan)
            limits2.append(np.nan)

            def model(model_norm_f, model_norm_r, model_ramp_a1,
                      model_ldc1, model_ldc2, model_ldc3, model_ldc4,
                      model_rp_over_rs, model_fp_over_fs, model_period, model_sma_over_rs, model_eccentricity,
                      model_inclination, model_periastron, model_mid_time, independent=False):

                normalisation = np.where(data_scan > 0, model_norm_f, model_norm_r)
                detrend1 = 1.0 - model_ramp_a1 * (data_time - model_mid_time)
                detrend2 = model_white

                if observation_type == 'transit':
                    signal = plc.transit(method.value, [model_ldc1, model_ldc2, model_ldc3, model_ldc4],
                                         model_rp_over_rs, model_period, model_sma_over_rs, model_eccentricity,
                                         model_inclination, model_periastron, model_mid_time, data_time_hr)
                else:
                    signal = plc.eclipse(model_fp_over_fs, model_rp_over_rs, model_period,
                                         model_sma_over_rs, model_eccentricity, model_inclination,
                                         model_periastron, model_mid_time + period.value / 2, data_time_hr)

                signal = np.mean(np.reshape(signal, (sub_exposures, len(data_time))), 0)

                if not independent:
                    return normalisation * detrend1 * detrend2 * signal
                else:
                    return normalisation * detrend1 * detrend2, signal

            curve_fit_fitted_parameters_indices = np.where(~np.isnan(np.array(limits1) * np.array(limits2)))[0]

            def curve_fit_model(xxx, *curve_fit_fitted_parameters):
                if xxx:
                    curve_fit_parameters = [fff for fff in initial]
                    for ii in range(len(curve_fit_fitted_parameters_indices)):
                        curve_fit_parameters[curve_fit_fitted_parameters_indices[ii]] = curve_fit_fitted_parameters[ii]

                    return model(*curve_fit_parameters)

            popt, pcov = curve_fit(curve_fit_model, 1, data_bin,
                                   p0=np.array(initial)[curve_fit_fitted_parameters_indices])

            data_bin_error *= (np.std(data_bin - curve_fit_model(1, *popt)) /
                               np.median(data_bin_error))

            mcmc = plc.EmceeFitting([data_bin, data_bin_error], model, np.array(initial), np.array(limits1),
                                    np.array(limits2), spectral_mcmc_walkers.value, spectral_mcmc_iterations.value,
                                    spectral_mcmc_burned_iterations.value, names, print_names, counter=False)

            mcmc.run_mcmc()

            spectral_fit = mcmc.results

            lower_wavelength.from_dictionary(bins_dictionaries[j])
            upper_wavelength.from_dictionary(bins_dictionaries[j])

            spectral_fit['wavelength'] = {'lambda1': lower_wavelength.value,
                                          'lambda2': upper_wavelength.value,
                                          'lambda_mean': 0.5 * (lower_wavelength.value + upper_wavelength.value),
                                          'lambda_width': upper_wavelength.value - lower_wavelength.value
                                          }

            spectral_fit['limb_darkening'] = {'method': method.value}

            spectral_fit['exposure'] = {'exp_time': exposure_time.value * (24.0 * 60 * 60),
                                        'model_resolution': sub_exposure_time * (24.0 * 60 * 60)
                                        }

            spectral_fit['input_time_series'] = {'hjd': data_time,
                                                 'scan': data_scan,
                                                 # 'white_raw_lc': data_white,
                                                 # 'white_raw_lc_error': data_white_error,
                                                 # 'white_model': model_white,
                                                 'raw_lc': data_bin,
                                                 'raw_lc_error': data_bin_error,
                                                 'raw_lc_ph_error': data_bin_ph_error
                                                 # 'relative_lc': data_relative,
                                                 # 'relative_lc_error': data_relative_error
                                                 }

            del spectral_fit['input_series']

            model_systematics, model_curve = model(*mcmc.results['parameters_final'], independent=True)
            spectral_fit['output_time_series'] = {'phase': ((data_time - spectral_fit['parameters']['t_0']['value']) /
                                                            spectral_fit['parameters']['P']['value']),
                                                  'full_model': model_systematics * model_curve,
                                                  'full_residuals': data_bin - model_systematics * model_curve,
                                                  'systematics': model_systematics,
                                                  'detrended_lc': data_bin / model_systematics,
                                                  'transit': model_curve,
                                                  'residuals': data_bin / model_systematics - model_curve
                                                  #
                                                  # 'systematics': model_systematics,
                                                  # 'detrended_lc': data_relative / model_systematics,
                                                  # 'transit': model_curve,
                                                  # 'residuals': data_relative / model_systematics - model_curve
                                                  }

            del spectral_fit['output_series']

            chis = np.round(
                np.sum((spectral_fit['output_time_series']['full_residuals'] ** 2) / (
                        spectral_fit['input_time_series']['raw_lc_error'] ** 2)), 1)
            rchis = np.round(
                np.sum((spectral_fit['output_time_series']['full_residuals'] ** 2) / (
                        spectral_fit['input_time_series']['raw_lc_error'] ** 2)) / (
                        len(data_time) - len(spectral_fit['statistics']['corr_matrix'][0])), 2)
            std = int(np.round(np.std(spectral_fit['output_time_series']['residuals']) * 1000000, 0))
            rstd = np.round(np.std(spectral_fit['output_time_series']['full_residuals']) / (
                np.mean([np.std(np.random.normal(0, spectral_fit['input_time_series']['raw_lc_ph_error']))
                         for ff in range(1000)])), 2)
            autocor = np.round(np.max(np.abs(spectral_fit['statistics']['res_autocorr'][1:3])), 2)
            spectral_fit['res_statistics'] = {'chisqr': chis,
                                              'rchisqr': rchis,
                                              'stdppm': std,
                                              'rstd': rstd,
                                              'autocor': autocor
                                              }

            fitting_results['lightcurves']['bin_{0}'.format(str(j + 1).zfill(2))] = spectral_fit

            counter.update()

        if observation_type == 'transit':
            lambda_mean = []
            rprs = []
            rprs_error = []
            lambda_width = []
            for j in range(len(bins_dictionaries)):
                bin_name = 'bin_{0}'.format(str(j + 1).zfill(2))
                lambda_mean.append(float(fitting_results['lightcurves'][bin_name]['wavelength']['lambda_mean'])
                                   / 10000.0)
                rprs.append(fitting_results['lightcurves'][bin_name]['parameters']['rp']['value'])
                rprs_error.append(np.mean([fitting_results['lightcurves'][bin_name]['parameters']['rp']['m_error'],
                                           fitting_results['lightcurves'][bin_name]['parameters']['rp']['p_error']]))
                lambda_width.append(float(fitting_results['lightcurves'][bin_name]['wavelength']['lambda_width'])
                                    / 10000.0)

            fitting_results['spectrum'] = {'wavelength': np.array(lambda_mean),
                                           'depth': np.array(rprs) ** 2,
                                           'error': 2 * np.array(rprs) * np.array(rprs_error),
                                           'width': np.array(lambda_width),
                                           }
        else:
            lambda_mean = []
            fpfs = []
            fpfs_error = []
            lambda_width = []
            for j in range(len(bins_dictionaries)):
                bin_name = 'bin_{0}'.format(str(j + 1).zfill(2))
                lambda_mean.append(float(fitting_results['lightcurves'][bin_name]['wavelength']['lambda_mean'])
                                   / 10000.0)
                fpfs.append(fitting_results['lightcurves'][bin_name]['parameters']['fp']['value'])
                fpfs_error.append(np.mean([fitting_results['lightcurves'][bin_name]['parameters']['fp']['m_error'],
                                           fitting_results['lightcurves'][bin_name]['parameters']['fp']['p_error']]))
                lambda_width.append(float(fitting_results['lightcurves'][bin_name]['wavelength']['lambda_width'])
                                    / 10000.0)
            fitting_results['spectrum'] = {'wavelength': np.array(lambda_mean),
                                           'fpfs': np.array(fpfs),
                                           'error': np.array(fpfs_error),
                                           'width': np.array(lambda_width),
                                           }

    return fitting_results


def plot_fitting(dictionary, directory):

    forward_colour = 'k'
    reverse_colour = 'r'

    def correlation(x, y):
            n = len(x)
            mx = np.mean(x)
            sx = np.std(x)
            my = np.mean(y)
            sy = np.std(y)
            return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)

    def od_distribution(data):

        data = np.array(data)
        xstep = np.sqrt(np.median((data - np.median(data)) ** 2)) / 5.0
        xmin = min(data)
        xmax = max(data)
        x_size = round((xmax - xmin) / xstep) + 1

        distrx = xmin + np.arange(x_size) * xstep
        data = np.int_((data - xmin) / xstep)
        distr = np.bincount(data)
        distr = np.insert(distr, len(distr), np.zeros(int(x_size) - len(distr)))

        plt.step(distrx, distr, c='k')

    def td_distribution(datax, datay):

        datax = np.array(datax)
        median = np.median(datax)
        med = np.sqrt(np.median((datax - median) ** 2))
        xstep = med / 5.0
        xmin = min(datax)
        xmax = max(datax)
        x_size = int(round((xmax - xmin) / xstep)) + 1
        datax = np.int_((datax - xmin) / xstep)
        datay = np.array(datay)
        median = np.median(datay)
        med = np.sqrt(np.median((datay - median) ** 2))
        ystep = med / 5.0
        ymin = min(datay)
        ymax = max(datay)
        y_size = int(round((ymax - ymin) / ystep)) + 1
        datay = np.int_((datay - ymin) / ystep)

        yx_size = x_size * y_size
        yx = datay * x_size + datax

        yx = np.bincount(yx)
        yx = np.insert(yx, len(yx), np.zeros(yx_size - len(yx)))

        xx, yy = np.meshgrid(xmin + np.arange(x_size) * xstep, ymin + np.arange(y_size) * ystep)

        final = np.reshape(yx, (y_size, x_size))
        plt.imshow(np.where(final > 0, np.log(np.where(final > 0, final, 1)), 0),
                   extent=(np.min(xx), np.max(xx), np.min(yy), np.max(yy)),
                   cmap=plt.cm.Greys, origin='lower', aspect='auto')

    def plot_correlations(light_curve_dic, bin_to_plot, export_file):

        names = []
        results = []
        errors1 = []
        errors2 = []
        errors = []
        traces = []

        for i in light_curve_dic[bin_to_plot]['statistics']['corr_variables'].split(','):
            names.append(light_curve_dic[bin_to_plot]['parameters'][i]['print_name'])
            results.append(light_curve_dic[bin_to_plot]['parameters'][i]['value'])
            errors1.append(light_curve_dic[bin_to_plot]['parameters'][i]['m_error'])
            errors2.append(light_curve_dic[bin_to_plot]['parameters'][i]['p_error'])
            errors.append(0.5 * (light_curve_dic[bin_to_plot]['parameters'][i]['m_error'] +
                                 light_curve_dic[bin_to_plot]['parameters'][i]['p_error']))
            traces.append(light_curve_dic[bin_to_plot]['parameters'][i]['trace'])

        all_var = len(traces)
        fig = plt.figure(figsize=(2.5 * all_var, 2.5 * all_var))
        fig.set_tight_layout(False)
        cmap = matplotlib.cm.get_cmap('brg')

        for var in range(len(names)):

            try:
                plt.subplot(all_var, all_var, all_var * var + var + 1, facecolor='w')
            except AttributeError:
                plt.subplot(all_var, all_var, all_var * var + var + 1, axisbg='w')

            od_distribution(traces[var])

            plt.axvline(results[var], c='k')
            plt.axvline(results[var] - errors1[var], c='k', ls='--', lw=0.5)
            plt.axvline(results[var] + errors2[var], c='k', ls='--', lw=0.5)

            plt.xticks(plt.xticks()[0], np.ones_like(plt.yticks()[0]))
            plt.yticks(plt.yticks()[0], np.ones_like(plt.yticks()[0]))
            plt.tick_params(left=False, right=False, top=False, bottom=False, labelbottom=False, labelleft=False)

            try:
                digit1 = abs(int(np.log10(errors1[var]))) + 1
            except OverflowError:
                digit1 = 3
            try:
                digit2 = abs(int(np.log10(errors2[var]))) + 1
            except OverflowError:
                digit2 = 3
            try:
                done1 = 1 // int(errors1[var] * (10 ** digit1))
            except ZeroDivisionError:
                done1 = 0
            try:
                done2 = 1 // int(errors2[var] * (10 ** digit2))
            except ZeroDivisionError:
                done2 = 0

            if errors1[var] > 1. and errors2[var] > 1.:
                digit1, digit2, done1, done2 = 0, 0, 0, 0

            width = max(digit1 + done1, digit2 + done2)

            plt.xlabel('{0}{1}{2}{3}{4}{5}{6}{7}{8}'.format(
                r'${0}$'.format(names[var]), '\n',
                r'${0:.{width}f}$'.format(round(results[var], width), width=width), '\n',
                r'$-$', r'${0:.{width}f}$'.format(round(errors1[var], width), width=width), '\n',
                r'$+$', r'${0:.{width}f}$'.format(round(errors2[var], width), width=width)),
                fontsize=15)

            plt.xlim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
            plt.ylim(0, plt.ylim()[1])

            for j in range(var + 1, all_var):

                plt.subplot(all_var, all_var, all_var * var + 1 + j)
                td_distribution(traces[j], traces[var])

                plt.yticks(plt.yticks()[0], np.arange(len(plt.yticks()[0])))
                plt.xticks(plt.xticks()[0], np.arange(len(plt.xticks()[0])))
                plt.tick_params(bottom=False, left=False, right=False, top=False, labelbottom=False,
                                labelleft=False, labelright=False, labeltop=False)

                plt.xlim(results[j] - 6 * errors[j], results[j] + 6 * errors[j])
                plt.ylim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
                text_x = plt.xlim()[1] - 0.05 * (plt.xlim()[1] - plt.xlim()[0])
                text_y = plt.ylim()[1] - 0.05 * (plt.ylim()[1] - plt.ylim()[0])
                plt.text(text_x, text_y, r'{0}{1}{2}'.format('$', str(correlation(traces[j], traces[var])), '$'),
                         color=cmap(abs(correlation(traces[j], traces[var])) / 2.), fontsize=15, ha='right', va='top')

        plt.subplots_adjust(hspace=0, wspace=0)
        functions.save_figure(directory, name=export_file)
        plt.close('all')

    def plot_fitting_i(light_curve_dic, bin_to_plot, export_file):

        forward_colour_i = 'k'
        reverse_colour_i = 'r'

        data_phase = light_curve_dic[bin_to_plot]['output_time_series']['phase']
        data_full = light_curve_dic[bin_to_plot]['input_time_series']['raw_lc']
        data_curve = light_curve_dic[bin_to_plot]['output_time_series']['detrended_lc']
        model_curve = light_curve_dic[bin_to_plot]['output_time_series']['transit']
        data_residuals = light_curve_dic[bin_to_plot]['output_time_series']['residuals']

        data_scan = np.array(light_curve_dic[bin_to_plot]['input_time_series']['scan'])
        fr = np.where(data_scan > 0)
        rv = np.where(data_scan < 0)

        fig = plt.figure(figsize=(10, 10))
        fig.set_tight_layout(False)

        plt.subplot(4, 1, 1)

        if len(fr[0]) > 0:
            plt.plot(data_phase[fr], data_full[fr] / np.mean(data_full[fr]), 'o', c=forward_colour_i,
                     mec=forward_colour_i, ms=3)
        if len(rv[0]) > 0:
            plt.plot(data_phase[rv], data_full[rv] / np.mean(data_full[rv]), 'o', c=reverse_colour_i,
                     mec=reverse_colour_i, ms=3)
        x_max = max(np.abs(plt.xlim()))

        plt.ylabel('{0}{1}{2}'.format('$\mathrm{raw}$', '\n', '$\mathrm{norm.} \, \mathrm{flux}$', fontsize=15))

        plt.xlim(-x_max, x_max)

        functions.adjust_ticks()
        plt.tick_params(labelbottom=False)

        plt.subplot(4, 1, 2)

        plt.plot(data_phase[fr], data_curve[fr], 'o', c=forward_colour_i, mec=forward_colour_i, ms=3)
        plt.plot(data_phase[rv], data_curve[rv], 'o', c=reverse_colour_i, mec=reverse_colour_i, ms=3)
        plt.plot(data_phase, model_curve, c='k', ls='--')

        plt.ylabel('{0}{1}{2}'.format('$\mathrm{de-trended}$', '\n', r'$\mathrm{norm.} \, \mathrm{flux}$', fontsize=15))

        plt.xlim(-x_max, x_max)

        functions.adjust_ticks()
        plt.tick_params(labelbottom=False)

        plt.subplot(4, 1, 3)

        plt.plot(data_phase[fr], (10 ** 6) * data_residuals[fr], 'o', c=forward_colour_i, mec=forward_colour_i, ms=3)
        plt.plot(data_phase[rv], (10 ** 6) * data_residuals[rv], 'o', c=reverse_colour_i, mec=reverse_colour_i, ms=3)
        plt.axhline(0, c='k', ls='--')

        plt.ylabel(r'$\mathrm{residuals} \, \mathrm{(ppm)}$', fontsize=15)
        plt.xlabel(r'$\mathrm{phase}$', fontsize=15)

        plt.xlim(-x_max, x_max)
        functions.adjust_ticks()

        plt.subplot(6, 1, 6)

        res_autocorr = np.array(light_curve_dic[bin_to_plot]['statistics']['res_autocorr'])
        plt.bar(np.arange(len(res_autocorr)), res_autocorr, width=0.5, color='k')

        plt.ylabel(r'$\mathrm{res.} \ \mathrm{ACF}$', fontsize=15)

        plt.ylim(-1.0, 1.0)
        functions.adjust_ticks()

        plt.subplots_adjust(hspace=0.0)
        functions.save_figure(directory, name=export_file)
        plt.close('all')

    def plot_fitting_all(light_curve_dic, export_file):

        dy = 0
        fig = plt.figure(figsize=(12, 10))
        fig.set_tight_layout(False)
        ax2 = plt.subplot(1, 3, 1)
        ax3 = plt.subplot(1, 3, 2)
        ax4 = plt.subplot(1, 3, 3)

        cmap = matplotlib.cm.get_cmap('rainbow')

        ii = 0
        ax4xlim = 0

        for ii in range(len(light_curve_dic) - 1):

            bin_name = 'bin_{0}'.format(str(ii + 1).zfill(2))

            bin_lc = light_curve_dic[bin_name]

            wavelength_mean = bin_lc['wavelength']['lambda_mean']
            wavelength_mean = 1.0 - (17000 - wavelength_mean) / (17000 - 7000)
            phase = np.array(bin_lc['output_time_series']['phase'])
            detrended_lc = np.array(bin_lc['output_time_series']['detrended_lc'])
            model_curve = np.array(bin_lc['output_time_series']['transit'])
            residuals = np.array(bin_lc['output_time_series']['residuals'])
            res_autocorr = np.array(bin_lc['statistics']['res_autocorr'])

            chis = bin_lc['res_statistics']['chisqr']
            rchis = bin_lc['res_statistics']['rchisqr']
            std = bin_lc['res_statistics']['stdppm']
            rstd = bin_lc['res_statistics']['rstd']
            autocor = bin_lc['res_statistics']['autocor']

            if dy == 0:
                dy = -8. * np.std(bin_lc['output_time_series']['residuals'])

            ax2.plot(phase, model_curve + ii * dy, c='k', ls='--', lw=0.5)
            ax2.plot(phase, detrended_lc + ii * dy, 'o', c=cmap(wavelength_mean), mew=0.1, ms=3)

            ax3.axhline(1 + ii * dy, c='k', ls='--', lw=0.5)
            ax3.plot(phase, 1 + residuals + ii * dy, 'o', c=cmap(wavelength_mean), mew=0.1, ms=3)

            ax4.axhline(1 + ii * dy, c='k', ls='--', lw=0.5)
            ax4.bar(np.arange(len(res_autocorr)), res_autocorr * (0.7 * -dy), bottom=1 + ii * dy,
                    linewidth=0.3, color=cmap(wavelength_mean))
            ax4xlim = ax4.get_xlim()[1]
            ax4.text(ax4xlim + 5, 1 + ii * dy, r'$\sigma$={0}, $\overline{1}$={2}, ${3}^2$={4}, '
                                                     r'$\overline{3}^2$={5}, $R^2$={6}'.format(
                std, '{\sigma}', rstd, '{\chi}', chis, rchis, autocor), va='center', ha='left')

        bin_lc = light_curve_dic['white']
        ii += 1
        phase = np.array(bin_lc['output_time_series']['phase'])
        detrended_lc = np.array(bin_lc['output_time_series']['detrended_lc'])
        model_curve = np.array(bin_lc['output_time_series']['transit'])
        residuals = np.array(bin_lc['output_time_series']['residuals'])
        res_autocorr = np.array(bin_lc['statistics']['res_autocorr'])

        chis = bin_lc['res_statistics']['chisqr']
        rchis = bin_lc['res_statistics']['rchisqr']
        std = bin_lc['res_statistics']['stdppm']
        rstd = bin_lc['res_statistics']['rstd']
        autocor = bin_lc['res_statistics']['autocor']

        ax2.plot(phase, model_curve + ii * dy, c='k', ls='--', lw=0.5)
        ax2.plot(phase, detrended_lc + ii * dy, 'o', c='k', mew=0.1, ms=3)

        x_max = max(np.abs(ax2.get_xlim()))
        y_min, y_max = ax2.get_ylim()

        ax3.axhline(1 + ii * dy, c='k', ls='--', lw=0.5)
        ax3.plot(phase, 1 + residuals + ii * dy, 'o', c='k', mew=0.1, ms=3)

        ax4.axhline(1 + ii * dy, c='k', ls='--', lw=0.5)
        ax4.bar(np.arange(len(res_autocorr)), res_autocorr * (0.7 * -dy), bottom=1 + ii * dy, linewidth=0.3, color='k')
        ax4.text(ax4.get_xlim()[1], 1 + ii * dy, r'$\sigma$={0}, $\overline{1}$={2}, ${3}^2$={4}, '
                                                 r'$\overline{3}^2$={5}, $R^2$={6}'.format(
            std, '{\sigma}', rstd, '{\chi}', chis, rchis, autocor), va='center', ha='left')

        ax2.set_ylabel(r'$\mathrm{norm.} \, \mathrm{flux}$', fontsize=15)
        ax2.set_title(r'$\mathrm{de-trended}$', fontsize=15)
        ax2.set_xlabel(r'$\mathrm{phase}$', fontsize=15)

        ax2.set_xlim(-x_max, x_max)
        ax2.set_ylim(y_min, y_max)
        ax2.set_xticklabels(ax2.get_xticks(), rotation=45, ha='right')
        functions.adjust_ticks_ax(ax2)

        ax3.set_title(r'$\mathrm{residuals}$', fontsize=15)

        ax3.set_xlim(-x_max, x_max)
        ax3.set_ylim(y_min, y_max)
        ax3.set_xticklabels(ax3.get_xticks(), rotation=45, ha='right')
        ax3.tick_params(labelleft=False)
        functions.adjust_ticks_ax(ax3)

        ax4.set_title(r'$\mathrm{res.} \ \mathrm{ACF}$', fontsize=15)
        ax4.set_xlabel(r'$\#$', fontsize=15)

        ax4.set_ylim(y_min, y_max)
        ax4.set_xlim(-5, ax4xlim)
        ax4.set_xticklabels(ax4.get_xticks(), rotation=45, ha='right')
        ax4.tick_params(labelleft=False)
        functions.adjust_ticks_ax(ax4)
        ax4.set_xlim(-5, ax4xlim)

        plt.subplots_adjust(wspace=0)
        functions.save_figure(directory, name=export_file)
        plt.close('all')

    def plot_spectral_results(lightcurve, export_file):

        wavelength_mean = []
        wavelength_width = []
        a1 = []
        a2 = []
        a3 = []
        a4 = []
        rp = []
        rp_er = []
        n_l_for = []
        n_l_for_er = []
        n_l_rev = []
        n_l_rev_er = []
        r_a1 = []
        r_a1_er = []

        eclipse = False

        for ii in range(len(lightcurve) - 1):
            bin_name = 'bin_{0}'.format(str(ii + 1).zfill(2))
            wavelength_mean.append(lightcurve[bin_name]['wavelength']['lambda_mean'])
            wavelength_width.append(lightcurve[bin_name]['wavelength']['lambda_width'])
            a1.append(lightcurve[bin_name]['parameters']['ldc_1']['value'])
            a2.append(lightcurve[bin_name]['parameters']['ldc_2']['value'])
            a3.append(lightcurve[bin_name]['parameters']['ldc_3']['value'])
            a4.append(lightcurve[bin_name]['parameters']['ldc_4']['value'])
            if lightcurve[bin_name]['parameters']['rp']['m_error'] is None:
                rp.append(lightcurve[bin_name]['parameters']['fp']['value'])
                rp_er.append(max(lightcurve[bin_name]['parameters']['fp']['m_error'],
                             lightcurve[bin_name]['parameters']['fp']['p_error']))
                eclipse = True
            else:
                rp.append(lightcurve[bin_name]['parameters']['rp']['value'])
                rp_er.append(max(lightcurve[bin_name]['parameters']['rp']['m_error'],
                                 lightcurve[bin_name]['parameters']['rp']['p_error']))
                eclipse = False
            if lightcurve[bin_name]['parameters']['n_l_for']['m_error']:
                n_l_for.append(lightcurve[bin_name]['parameters']['n_l_for']['value'])
                n_l_for_er.append(max(lightcurve[bin_name]['parameters']['n_l_for']['m_error'],
                                      lightcurve[bin_name]['parameters']['n_l_for']['p_error']))
            if lightcurve[bin_name]['parameters']['n_l_rev']['m_error']:
                n_l_rev.append(lightcurve[bin_name]['parameters']['n_l_rev']['value'])
                n_l_rev_er.append(max(lightcurve[bin_name]['parameters']['n_l_rev']['m_error'],
                                      lightcurve[bin_name]['parameters']['n_l_rev']['p_error']))

            r_a1.append(lightcurve[bin_name]['parameters']['r_a1']['value'])
            r_a1_er.append(max(lightcurve[bin_name]['parameters']['r_a1']['m_error'],
                               lightcurve[bin_name]['parameters']['r_a1']['p_error']))

        wavelength_mean = np.array(wavelength_mean) / 10000
        wavelength_width = np.array(wavelength_width) / 10000
        a1 = np.array(a1)
        a2 = np.array(a2)
        a3 = np.array(a3)
        a4 = np.array(a4)
        rp = np.array(rp)
        rp_er = np.array(rp_er)
        n_l_for = np.array(n_l_for)
        n_l_for_er = np.array(n_l_for_er)
        n_l_rev = np.array(n_l_rev)
        n_l_rev_er = np.array(n_l_rev_er)

        # plot inputs and results

        if len(n_l_rev) > 0 and len(n_l_for) > 0:
            n_plots = 6
        else:
            n_plots = 5

        fig = plt.figure(figsize=(10, n_plots * 2.5))
        fig.set_tight_layout(False)
        ax1 = plt.subplot(n_plots, 1, 1)
        ax1.cla()
        ax2 = plt.subplot(n_plots, 1, 2)
        ax2.cla()
        ax3 = plt.subplot(n_plots, 1, 3)
        ax3.cla()
        ax4 = plt.subplot(n_plots, 1, 4)
        ax4.cla()
        ax5 = plt.subplot(n_plots, 1, 5)
        ax5.cla()

        # plot 1

        ax1.errorbar(wavelength_mean, wavelength_width, xerr=wavelength_width / 2.,
                     fmt='o', color='k', ms=n_plots, mec='None')
        ax1.set_ylabel(r'$\mathrm{bin} \ \mathrm{width} \, (\mu \mathrm{m})$', fontsize=15)
        ax1.tick_params(labelbottom=False)
        functions.adjust_ticks_ax(ax1)

        # plot 2

        ax2.plot(wavelength_mean, a1, 'bo', ms=5, mec='None', label=r'$\mathrm{ldc_1}$')
        ax2.plot(wavelength_mean, a2, 'ro', ms=5, mec='None', label=r'$\mathrm{ldc_2}$')
        ax2.plot(wavelength_mean, a3, 'go', ms=5, mec='None', label=r'$\mathrm{ldc_3}$')
        ax2.plot(wavelength_mean, a4, 'ko', ms=5, mec='None', label=r'$\mathrm{ldc_4}$')
        ax2.legend()
        ax2.set_ylabel('{0}{1}{2}'.format(r'$\mathrm{limb}$', '\n', r'$\mathrm{darkening}$', fontsize=15))
        ax2.tick_params(labelbottom=False)
        functions.adjust_ticks_ax(ax2)

        # plot 3

        ax3.errorbar(wavelength_mean, rp, rp_er, color='k', fmt='-o', mec='None')

        if eclipse:
            ax3.set_ylabel(r'$F_\mathrm{p} / F_*}$', fontsize=15)
        else:
            ax3.set_ylabel(r'$R_\mathrm{p} / R_*}$', fontsize=15)

        dy = 0.6 * (np.max(rp) - np.min(rp)) + np.max(rp_er)
        ax3.set_ylim((np.max(rp) + np.min(rp)) / 2 - dy, (np.max(rp) + np.min(rp)) / 2 + dy)
        ax3.tick_params(labelbottom=False)
        functions.adjust_ticks_ax(ax3)

        # plot 4

        next_ax = ax5

        if len(n_l_for) > 0 and n_plots == 5:

            ax4.errorbar(wavelength_mean, n_l_for, n_l_for_er, color='k', fmt='-o', mec='None')
            ax4.set_ylabel(r'$n_\lambda ^\mathrm{for}$', fontsize=15)
            ax4.tick_params(labelbottom=False)
            functions.adjust_ticks_ax(ax4)

        elif len(n_l_rev) > 0 and n_plots == 5:

            ax4.errorbar(wavelength_mean, n_l_rev, n_l_rev_er, color='k', fmt='-o', mec='None')
            ax4.set_ylabel(r'$n_\lambda ^\mathrm{for}$', fontsize=15)
            ax4.tick_params(labelbottom=False)
            functions.adjust_ticks_ax(ax4)

        elif n_plots == 6:

            ax4.errorbar(wavelength_mean, n_l_for, n_l_for_er, color='k', fmt='-o', mec='None')
            ax4.set_ylabel(r'$n_\lambda ^\mathrm{for}$', fontsize=15)
            ax4.tick_params(labelbottom=False)
            functions.adjust_ticks_ax(ax4)

            ax5.errorbar(wavelength_mean, n_l_rev, n_l_rev_er, color='k', fmt='-o', mec='None')
            ax5.set_ylabel(r'$n_\lambda ^\mathrm{for}$', fontsize=15)
            ax5.tick_params(labelbottom=False)
            functions.adjust_ticks_ax(ax5)

            next_ax = plt.subplot(n_plots, 1, 6)

        # plot 5

        next_ax.errorbar(wavelength_mean, r_a1, r_a1_er, color='k', fmt='-o', mec='None')

        next_ax.set_ylabel(r'$\mathrm{ramp} \ \mathrm{slope}$', fontsize=15)
        next_ax.set_xlabel(r'$\lambda \, (\mu \mathrm{m})$', fontsize=15)

        functions.adjust_ticks_ax(next_ax)

        plt.subplots_adjust(hspace=0)
        functions.save_figure(directory, name=export_file)
        plt.close('all')

    def plot_diagnostics(lightcurve, export_file):

        fig = plt.figure(10, figsize=(7, 10))
        fig.set_tight_layout(False)

        hjd_time = lightcurve['white']['input_time_series']['hjd']
        flux = lightcurve['white']['input_time_series']['raw_lc']
        yshift = lightcurve['white']['input_time_series']['y_shift']
        yshift_err = lightcurve['white']['input_time_series']['y_shift_error']
        xshift = lightcurve['white']['input_time_series']['x_shift']
        xshift_err = lightcurve['white']['input_time_series']['x_shift_error']
        scan = lightcurve['white']['input_time_series']['scan']
        ssky = lightcurve['white']['input_time_series']['sky']

        scan = np.array(scan)
        forward = np.where(scan > 0)
        reverse = np.where(scan < 0)
        if len(forward[0]) > 0:
            for_flag = True
        else:
            for_flag = False
        if len(reverse[0]) > 0:
            rev_flag = True
        else:
            rev_flag = False

        hjd_time = np.array(hjd_time) - hjd_time[0]
        xshift = np.array(xshift)
        if for_flag:
            xshift[forward] = xshift[forward] - xshift[forward][0]
        if rev_flag:
            xshift[reverse] = xshift[reverse] - xshift[reverse][0]
        xshift_err = np.array(xshift_err)
        yshift = np.array(yshift)
        if for_flag:
            yshift[forward] = yshift[forward] - yshift[forward][0]
        if rev_flag:
            yshift[reverse] = yshift[reverse] - yshift[reverse][0]
        yshift_err = np.array(yshift_err)
        ssky = np.array(ssky)

        plt.subplot(4, 1, 1)
        plt.cla()
        plt.errorbar(hjd_time[forward], yshift[forward], yshift_err[forward],
                     fmt='o', c=forward_colour, mec=forward_colour, ms=3)
        plt.errorbar(hjd_time[reverse], yshift[reverse], yshift_err[reverse],
                     fmt='o', c=reverse_colour, mec=reverse_colour, ms=3)

        plt.ylabel(r'$\Delta y_i \, \mathrm{(pix)}$', fontsize=15)

        functions.adjust_ticks()
        plt.tick_params(labelbottom=False)

        plt.subplot(4, 1, 2)
        plt.cla()
        plt.errorbar(hjd_time[forward], xshift[forward], xshift_err[forward],
                     fmt='o', c=forward_colour, mec=forward_colour, ms=3)
        plt.errorbar(hjd_time[reverse], xshift[reverse], xshift_err[reverse],
                     fmt='o', c=reverse_colour, mec=reverse_colour, ms=3)

        functions.adjust_ticks()
        plt.ylabel(r'$\Delta x_i \, \mathrm{(pix)}$', fontsize=15)
        plt.tick_params(labelbottom=False)

        plt.subplot(4, 1, 3)
        plt.cla()
        plt.plot(hjd_time[forward], ssky[forward],
                 'o', c=forward_colour, mec=forward_colour, ms=3)
        plt.plot(hjd_time[reverse], ssky[reverse],
                 'o', c=reverse_colour, mec=reverse_colour, ms=3)

        functions.adjust_ticks()
        plt.ylabel(r'$\mathrm{sky} \, \mathrm{ratio}$', fontsize=15)
        plt.tick_params(labelbottom=False)

        plt.subplot(4, 1, 4)
        plt.cla()
        plt.plot((np.array(hjd_time) - hjd_time[0])[forward], np.array(flux)[forward] / (10 ** 8),
                 'o', c=forward_colour, mec=forward_colour, ms=3)
        plt.plot((np.array(hjd_time) - hjd_time[0])[reverse], np.array(flux)[reverse] / (10 ** 8),
                 'o', c=reverse_colour, mec=reverse_colour, ms=3)

        functions.adjust_ticks()

        plt.ylabel(r'$\mathrm{e}^{-} \, (\times 10^8)$', fontsize=15)
        plt.xlabel(r'$\Delta t \, \mathrm{(days)}$', fontsize=15)

        plt.subplots_adjust(hspace=0)
        functions.save_figure(directory, name=export_file)
        plt.close('all')

    plot_diagnostics(dictionary['lightcurves'], 'diagnostics')
    plot_correlations(dictionary['lightcurves'], 'white', 'white_correlations')
    plot_correlations(dictionary['lightcurves'], 'white', 'white_correlations')
    plot_fitting_i(dictionary['lightcurves'], 'white', 'white_fitting')
    plot_fitting_all(dictionary['lightcurves'], 'all_fitting')

    # noinspection PyBroadException
    try:
        plot_spectral_results(dictionary['lightcurves'], 'spectral_results')
    except:
        pass

    try:
        plot_correlations(dictionary['lightcurves'], 'bin_10', 'bin_correlations')
    except KeyError:
        # noinspection PyBroadException
        try:
            plot_correlations(dictionary['lightcurves'], 'bin_01', 'bin_correlations')
        except:
            pass

    try:
        plot_fitting_i(dictionary['lightcurves'], 'bin_10', 'bin_fitting')
    except KeyError:
        # noinspection PyBroadException
        try:
            plot_fitting_i(dictionary['lightcurves'], 'bin_01', 'bin_fitting')
        except:
            pass
