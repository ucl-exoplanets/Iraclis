"""
images3_photometry.py

Includes all the functions that perform photometry processes.
All the functions take as input either an HDUList object or a DataSet object, as defined in the basics.py file and
return the input object and a dictionary that contains the extracted light-curves. In all cases, the default values for
the input parameters are the values in the respective pipeline.variables object. Note that the parameters for the
supporting functions do not have default values, as their purpose is to be used only in this particular file.

Functions included:
photometry:         ...
split_photometry:   ...

Supporting functions included:
get_flux_integral:  ...
get_flux_gauss:     ...

"""

from ._3objects import *


def get_flux_integral(fits, lower_wavelength, upper_wavelength,
                      aperture_lower_extend, aperture_upper_extend, sigma, plot=False):

    x_star = variables.x_star.custom_from_fits(fits).value
    y_star = variables.y_star.custom_from_fits(fits).value
    spectrum_direction = variables.spectrum_direction.custom_from_fits(fits).value
    scan_length = variables.scan_length.custom_from_fits(fits).value
    wdpt_constant_coefficient_1 = variables.wdpt_constant_coefficient_1.custom_from_fits(fits).value
    wdpt_constant_coefficient_2 = variables.wdpt_constant_coefficient_2.custom_from_fits(fits).value
    wdpt_constant_coefficient_3 = variables.wdpt_constant_coefficient_3.custom_from_fits(fits).value
    wdpt_slope_coefficient_1 = variables.wdpt_slope_coefficient_1.custom_from_fits(fits).value
    wdpt_slope_coefficient_2 = variables.wdpt_slope_coefficient_2.custom_from_fits(fits).value
    wdpt_slope_coefficient_3 = variables.wdpt_slope_coefficient_3.custom_from_fits(fits).value

    trace_at0 = calibrations.trace_at0.match(fits)
    trace_at1 = calibrations.trace_at1.match(fits)
    trace_at2 = calibrations.trace_at2.match(fits)
    trace_at3 = calibrations.trace_at3.match(fits)
    trace_at4 = calibrations.trace_at4.match(fits)
    trace_at5 = calibrations.trace_at5.match(fits)
    trace_bt0 = calibrations.trace_bt0.match(fits)
    trace_bt1 = calibrations.trace_bt1.match(fits)
    trace_bt2 = calibrations.trace_bt2.match(fits)

    def get_trace(dy):

        xx0 = x_star
        yy0 = y_star + dy
        sub = 507 - len(fits[1].data) / 2

        bt = trace_bt0 + trace_bt1 * xx0 + trace_bt2 * yy0
        at = (trace_at0 + trace_at1 * xx0 + trace_at2 * yy0 + trace_at3 * xx0 * xx0 +
              trace_at4 * xx0 * yy0 + trace_at5 * yy0 * yy0)
        return at, bt + yy0 - at * xx0 - sub + at * sub

    if spectrum_direction > 0:
        y0 = aperture_lower_extend
        y1 = scan_length + aperture_upper_extend
    else:
        y0 = - scan_length - aperture_upper_extend
        y1 = - aperture_lower_extend

    va1 = (wdpt_slope_coefficient_1 / (wdpt_slope_coefficient_2 + lower_wavelength) + wdpt_slope_coefficient_3)
    vb1 = (wdpt_constant_coefficient_1 / (wdpt_constant_coefficient_2 + lower_wavelength) + wdpt_constant_coefficient_3)
    va2 = (wdpt_slope_coefficient_1 / (wdpt_slope_coefficient_2 + upper_wavelength) + wdpt_slope_coefficient_3)
    vb2 = (wdpt_constant_coefficient_1 / (wdpt_constant_coefficient_2 + upper_wavelength) + wdpt_constant_coefficient_3)
    ha1, hb1 = get_trace(y0)
    ha2, hb2 = get_trace(y1)

    ha2 += sigma
    ha2 -= sigma

    if plot:
        xxx = np.arange((hb1 - vb1) / (va1 - ha1), (hb1 - vb2) / (va2 - ha1), 0.0001)
        plt.plot(xxx, ha1 * xxx + hb1, 'w-')
        xxx = np.arange((hb2 - vb1) / (va1 - ha2), (hb2 - vb2) / (va2 - ha2), 0.0001)
        plt.plot(xxx, ha2 * xxx + hb2, 'w-')
        xxx = np.arange((hb2 - vb1) / (va1 - ha2), (hb1 - vb1) / (va1 - ha1), 0.0001)
        plt.plot(xxx, va1 * xxx + vb1, 'w-')
        xxx = np.arange((hb2 - vb2) / (va2 - ha2), (hb1 - vb2) / (va2 - ha1), 0.0001)
        plt.plot(xxx, va2 * xxx + vb2, 'w-')

    fcr = np.full_like(fits[1].data, fits[1].data)
    fhm = np.roll(fcr, 1, axis=1)
    fhp = np.roll(fcr, -1, axis=1)
    fvm = np.roll(fcr, -1, axis=0)
    fvp = np.roll(fcr, 1, axis=0)
    x0, y0 = np.meshgrid(np.arange(len(fcr)), np.arange(len(fcr)))

    summ1 = (2.0 * fcr - fhm - fhp)
    summ2 = (4.0 * fcr - 4.0 * fhm)
    summ3 = (8.0 * fcr + 4.0 * fhm - 2.0 * fhp + 4.0 * fvm - 2.0 * fvp)
    summ4 = (4.0 * fcr - 4.0 * fvm)
    summ5 = (2.0 * fcr - fvm - fvp)

    summ6 = (4.0 * fcr - 4.0 * fhp)
    summ7 = (10.0 * fcr - 2.0 * fhm + 4.0 * fhp - fvm + fvp)
    summ8 = (20.0 * fcr - 4.0 * fhm + 8.0 * fhp)

    summ9 = (8.0 * fcr - 2.0 * fhm + 4.0 * fhp - 2.0 * fvm + 4.0 * fvp)
    summ10 = (4.0 * fcr - 4.0 * fvp)
    summ11 = (2.0 * fcr - fvm - fvp)

    # left edge
    a, b = va1, vb1
    x1 = (-b + y0) / a
    x2 = (1 - b + y0) / a

    formula = a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case1 = (
        + fcr - (1.0 / (24.0 * (a ** 3))) * (
            + summ1 * formula_4
            + a * summ2 * formula_3
            + (a ** 2) * (- summ3 * formula_2 - summ4 * formula_3 + summ5 * formula_4)
        )
    )

    formula = a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula

    case2 = (
        - (1.0 / (24.0 * (a ** 3))) * (
            + 4.0 * summ1 * (- 0.25 + formula - 1.5 * formula_2 + formula_3)
            + (a * 3.0) * summ6 * (-1.0 / 3 + formula - formula_2)
            + (a ** 2) * (summ7 - summ8 * formula)
        )
    )

    formula = - 1.0 + a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case3 = (
        - (1.0 / (24.0 * (a ** 3))) * (
            - summ1 * formula_4
            + a * summ6 * formula_3
            + (a ** 2) * (summ9 * formula_2 - summ10 * formula_3 - summ11 * formula_4)
        )
    )

    new_data = np.full_like(fits[1].data, fits[1].data)
    new_data = np.where((x1 > x0) & (x2 < x0), case1, new_data)
    new_data = np.where((x1 > x0) & (x0 + 1 > x1) & (x2 > x0) & (x0 + 1 > x2), case2, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 < x0 + 1), case3, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 > x0 + 1), 0, new_data)

    # right edge
    a, b = va2, vb2
    x1 = (-b + y0) / a
    x2 = (1 - b + y0) / a

    formula = a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case1 = (
        + (1.0 / (24.0 * (a ** 3))) * (
            + summ1 * formula_4
            + a * summ2 * formula_3
            + (a ** 2) * (- summ3 * formula_2 - summ4 * formula_3 + summ5 * formula_4)
        )
    )

    formula = a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula

    case2 = (
        fcr + (1.0 / (24.0 * (a ** 3))) * (
            + 4.0 * summ1 * (- 0.25 + formula - 1.5 * formula_2 + formula_3)
            + (a * 3.0) * summ6 * (-1.0 / 3 + formula - formula_2)
            + (a ** 2) * (summ7 - summ8 * formula)
        )
    )

    formula = - 1.0 + a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case3 = (
        fcr + (1.0 / (24.0 * (a ** 3))) * (
            - summ1 * formula_4
            + a * summ6 * formula_3
            + (a ** 2) * (summ9 * formula_2 - summ10 * formula_3 - summ11 * formula_4)
        )
    )

    new_data = np.where((x1 < x0) & (x2 < x0), 0, new_data)
    new_data = np.where((x1 > x0) & (x2 < x0), case1, new_data)
    new_data = np.where((x1 > x0) & (x0 + 1 > x1) & (x2 > x0) & (x0 + 1 > x2), case2, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 < x0 + 1), case3, new_data)

    # upper edge

    new_data = np.rot90(new_data)

    fcr = np.ones_like(new_data) * new_data
    fhm = np.roll(fcr, 1, axis=1)
    fhp = np.roll(fcr, -1, axis=1)
    fvm = np.roll(fcr, -1, axis=0)
    fvp = np.roll(fcr, 1, axis=0)
    x0, y0 = np.meshgrid(np.arange(len(fcr)), np.arange(len(fcr)))

    summ1 = (2.0 * fcr - fhm - fhp)
    summ2 = (4.0 * fcr - 4.0 * fhm)
    summ3 = (8.0 * fcr + 4.0 * fhm - 2.0 * fhp + 4.0 * fvm - 2.0 * fvp)
    summ4 = (4.0 * fcr - 4.0 * fvm)
    summ5 = (2.0 * fcr - fvm - fvp)

    summ6 = (4.0 * fcr - 4.0 * fhp)
    summ7 = (10.0 * fcr - 2.0 * fhm + 4.0 * fhp - fvm + fvp)
    summ8 = (20.0 * fcr - 4.0 * fhm + 8.0 * fhp)

    summ9 = (8.0 * fcr - 2.0 * fhm + 4.0 * fhp - 2.0 * fvm + 4.0 * fvp)
    summ10 = (4.0 * fcr - 4.0 * fvp)
    summ11 = (2.0 * fcr - fvm - fvp)

    a, b = ha2, hb2
    a, b = - 1.0 / a, len(fcr) + b / a
    x1 = (-b + y0) / a
    x2 = (1 - b + y0) / a

    formula = a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case1 = (
        + (1.0 / (24.0 * (a ** 3))) * (
            + summ1 * formula_4
            + a * summ2 * formula_3
            + (a ** 2) * (- summ3 * formula_2 - summ4 * formula_3 + summ5 * formula_4)
        )
    )

    formula = a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula

    case2 = (
        fcr + (1.0 / (24.0 * (a ** 3))) * (
            + 4.0 * summ1 * (- 0.25 + formula - 1.5 * formula_2 + formula_3)
            + (a * 3.0) * summ6 * (-1.0 / 3 + formula - formula_2)
            + (a ** 2) * (summ7 - summ8 * formula)
        )
    )

    formula = - 1.0 + a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case3 = (
        fcr + (1.0 / (24.0 * (a ** 3))) * (
            - summ1 * formula_4
            + a * summ6 * formula_3
            + (a ** 2) * (summ9 * formula_2 - summ10 * formula_3 - summ11 * formula_4)
        )
    )

    new_data = np.where((x1 < x0) & (x2 < x0), 0, new_data)
    new_data = np.where((x1 > x0) & (x2 < x0), case1, new_data)
    new_data = np.where((x1 > x0) & (x0 + 1 > x1) & (x2 > x0) & (x0 + 1 > x2), case2, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 < x0 + 1), case3, new_data)

    # lower edge

    a, b = ha1, hb1
    a, b = - 1.0 / a, len(fcr) + b / a
    x1 = (-b + y0) / a
    x2 = (1 - b + y0) / a

    formula = a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case1 = (
        + fcr - (1.0 / (24.0 * (a ** 3))) * (
            + summ1 * formula_4
            + a * summ2 * formula_3
            + (a ** 2) * (- summ3 * formula_2 - summ4 * formula_3 + summ5 * formula_4)
        )
    )

    formula = a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula

    case2 = (
        - (1.0 / (24.0 * (a ** 3))) * (
            + 4.0 * summ1 * (- 0.25 + formula - 1.5 * formula_2 + formula_3)
            + (a * 3.0) * summ6 * (-1.0 / 3 + formula - formula_2)
            + (a ** 2) * (summ7 - summ8 * formula)
        )
    )

    formula = - 1.0 + a + a * x0 + b - y0
    formula_2 = formula * formula
    formula_3 = formula_2 * formula
    formula_4 = formula_3 * formula

    case3 = (
        - (1.0 / (24.0 * (a ** 3))) * (
            - summ1 * formula_4
            + a * summ6 * formula_3
            + (a ** 2) * (summ9 * formula_2 - summ10 * formula_3 - summ11 * formula_4)
        )
    )

    new_data = np.where((x1 > x0) & (x2 < x0), case1, new_data)
    new_data = np.where((x1 > x0) & (x0 + 1 > x1) & (x2 > x0) & (x0 + 1 > x2), case2, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 < x0 + 1), case3, new_data)
    new_data = np.where((x1 > x0 + 1) & (x2 > x0 + 1), 0, new_data)

    new_data = np.rot90(new_data, 3)

    # error array

    xx = np.where(fits[1].data == 0, 1, fits[1].data)
    error = np.sqrt(new_data / xx) * fits[2].data

    flux = np.sum(new_data)
    error = np.sqrt(np.nansum(error * error))

    return flux, error


def get_flux_gauss(fits, lower_wavelength, upper_wavelength,
                   aperture_lower_extend, aperture_upper_extend, sigma, plot=False):

    spectrum_direction = variables.spectrum_direction.custom_from_fits(fits).value
    scan_length = variables.scan_length.custom_from_fits(fits).value
    scan_frame = variables.scan_frame.custom_from_fits(fits).value
    wavelength_frame = variables.wavelength_frame.custom_from_fits(fits).value

    if spectrum_direction > 0:
        y1 = min(aperture_lower_extend, aperture_upper_extend)
        y2 = scan_length + max(aperture_lower_extend, aperture_upper_extend)
    else:
        y1 = - scan_length - max(aperture_lower_extend, aperture_upper_extend)
        y2 = - min(aperture_lower_extend, aperture_upper_extend)

    science_frame = np.array(fits[plc.fits_sci(fits)[0]].data)
    error_frame = np.array(fits[plc.fits_err(fits)[0]].data)
    ph_error_frame = np.sqrt(np.abs(science_frame))

    scan_weight = (scipy.special.erf((scan_frame - y1) / ((sigma / 45.) * np.sqrt(2.0))) -
                   scipy.special.erf((scan_frame - y2) / ((sigma / 45.) * np.sqrt(2.0)))) / 2

    wavelength_weight = (scipy.special.erf((wavelength_frame - lower_wavelength) / (sigma * np.sqrt(2.0))) -
                         scipy.special.erf((wavelength_frame - upper_wavelength) / (sigma * np.sqrt(2.0)))) / 2

    weighted_science_frame = science_frame * scan_weight * wavelength_weight
    weighted_error_frame = error_frame * scan_weight * wavelength_weight
    weighted_ph_error_frame = ph_error_frame * scan_weight * wavelength_weight

    flux = np.sum(weighted_science_frame)
    error = np.sqrt(np.nansum(weighted_error_frame * weighted_error_frame))
    ph_error = np.sqrt(np.nansum(weighted_ph_error_frame * weighted_ph_error_frame))

    if plot:
        get_flux_integral(fits, lower_wavelength, upper_wavelength,
                          aperture_lower_extend, aperture_upper_extend, sigma, plot=True)

    return flux, error, ph_error


def photometry(input_data, white_lower_wavelength=None, white_upper_wavelength=None, bins_file=None,
               aperture_lower_extend=None, aperture_upper_extend=None, extraction_method=None,
               extraction_gauss_sigma=None, plot=False):

    # load pipeline and calibration variables to be used

    white_lower_wavelength = variables.white_lower_wavelength.custom(white_lower_wavelength)
    white_upper_wavelength = variables.white_upper_wavelength.custom(white_upper_wavelength)
    bins_file = variables.bins_file.custom(bins_file)
    aperture_lower_extend = variables.aperture_lower_extend.custom(aperture_lower_extend)
    aperture_upper_extend = variables.aperture_upper_extend.custom(aperture_upper_extend)
    extraction_method = variables.extraction_method.custom(extraction_method)
    extraction_gauss_sigma = variables.extraction_gauss_sigma.custom(extraction_gauss_sigma)

    ra_target = variables.ra_target.custom()
    dec_target = variables.dec_target.custom()
    subarray_size = variables.sub_array_size.custom()
    grism = variables.grism.custom()
    exposure_time = variables.exposure_time.custom()
    bins_number = variables.bins_number.custom()
    heliocentric_julian_date = variables.heliocentric_julian_date.custom()
    spectrum_direction = variables.spectrum_direction.custom()
    sky_background_level = variables.sky_background_level.custom()
    y_star = variables.y_star.custom()
    y_shift_error = variables.y_shift_error.custom()
    x_star = variables.x_star.custom()
    x_shift_error = variables.x_shift_error.custom()
    scan_length = variables.scan_length.custom()
    scan_length_error = variables.scan_length_error.custom()
    heliocentric_julian_date_array = variables.heliocentric_julian_date_array.custom()
    spectrum_direction_array = variables.spectrum_direction_array.custom()
    sky_background_level_array = variables.sky_background_level_array.custom()
    x_star_array = variables.x_star_array.custom()
    x_shift_error_array = variables.x_shift_error_array.custom()
    y_star_array = variables.y_star_array.custom()
    y_shift_error_array = variables.y_shift_error_array.custom()
    scan_length_array = variables.scan_length_array.custom()
    scan_length_error_array = variables.scan_length_error_array.custom()
    white_ldc1 = variables.white_ldc1.custom()
    white_ldc2 = variables.white_ldc2.custom()
    white_ldc3 = variables.white_ldc3.custom()
    white_ldc4 = variables.white_ldc4.custom()

    lower_wavelength = variables.lower_wavelength.custom()
    upper_wavelength = variables.upper_wavelength.custom()
    flux_array = variables.flux_array.custom()
    error_array = variables.error_array.custom()
    ph_error_array = variables.ph_error_array.custom()

    # set bins

    white_dictionary, bins_dictionaries = \
        variables.set_binning(input_data, white_lower_wavelength.value, white_upper_wavelength.value,
                              white_ldc1.value, white_ldc2.value, white_ldc3.value, white_ldc4.value,
                              bins_file.value)

    # select extraction method

    used_extraction_method = {'integral': get_flux_integral, 'gauss': get_flux_gauss}[extraction_method.value]

    # initiate counter

    counter = PipelineCounter('Photometry', len(input_data.spectroscopic_images))

    # iterate over the list of HDUList objects included in the input data

    light_curve = {}

    for fits in input_data.spectroscopic_images:

        try:
            ra_target.from_dictionary(light_curve)

        except KeyError:

            ra_target.from_fits(fits)
            ra_target.to_dictionary(light_curve)

            dec_target.from_fits(fits)
            dec_target.to_dictionary(light_curve)

            subarray_size.set(len(fits[1].data))
            subarray_size.to_dictionary(light_curve)

            grism.from_fits(fits)
            grism.to_dictionary(light_curve)

            exposure_time.from_fits(fits)
            exposure_time.to_dictionary(light_curve)

            aperture_lower_extend.to_dictionary(light_curve)

            aperture_upper_extend.to_dictionary(light_curve)

            extraction_method.to_dictionary(light_curve)

            extraction_gauss_sigma.to_dictionary(light_curve)

        heliocentric_julian_date.from_fits(fits)
        heliocentric_julian_date_array.set(
            np.append(heliocentric_julian_date_array.value, heliocentric_julian_date.value))

        spectrum_direction.from_fits(fits)
        spectrum_direction_array.set(np.append(spectrum_direction_array.value, spectrum_direction.value))

        sky_background_level.from_fits(fits, position=plc.fits_sci(fits)[0])
        sky_background_level_array.set(np.append(sky_background_level_array.value, sky_background_level.value))

        y_star.from_fits(fits)
        y_star_array.set(np.append(y_star_array.value, y_star.value))

        y_shift_error.from_fits(fits)
        y_shift_error_array.set(np.append(y_shift_error_array.value, y_shift_error.value))

        x_star.from_fits(fits)
        x_star_array.set(np.append(x_star_array.value, x_star.value))

        x_shift_error.from_fits(fits)
        x_shift_error_array.set(np.append(x_shift_error_array.value, x_shift_error.value))

        scan_length.from_fits(fits)
        scan_length_array.set(np.append(scan_length_array.value, scan_length.value))

        scan_length_error.from_fits(fits)
        scan_length_error_array.set(np.append(scan_length_error_array.value, scan_length_error.value))

        bins_number.set(len(bins_dictionaries))
        bins_number.to_dictionary(light_curve)

        for i in [white_dictionary] + bins_dictionaries:
            lower_wavelength.from_dictionary(i)
            upper_wavelength.from_dictionary(i)
            flux, error, ph_error = used_extraction_method(fits, lower_wavelength.value, upper_wavelength.value,
                                                           aperture_lower_extend.value, aperture_upper_extend.value,
                                                           extraction_gauss_sigma.value)
            flux_array.from_dictionary(i)
            flux_array.to_dictionary(i, value=np.append(flux_array.value, flux))
            error_array.from_dictionary(i)
            error_array.to_dictionary(i, value=np.append(error_array.value, error))
            ph_error_array.from_dictionary(i)
            ph_error_array.to_dictionary(i, value=np.append(ph_error_array.value, ph_error))

        counter.update()

        if plot:

            plt.figure(1)
            plt.imshow(fits[1].data, origin='lower', aspect='auto')
            plt.xlim(0, len(fits[1].data))
            plt.ylim(0, len(fits[1].data))

            lower_wavelength.from_dictionary(white_dictionary)
            upper_wavelength.from_dictionary(white_dictionary)

            used_extraction_method(fits, lower_wavelength.value, upper_wavelength.value,
                                   aperture_lower_extend.value, aperture_upper_extend.value,
                                   extraction_gauss_sigma.value, plot=True)

            plt.xlabel(r'$\mathrm{column \, (pix)}$', fontsize=20)
            plt.ylabel(r'$\mathrm{row \, (pix)}$', fontsize=20)

            plt.figure(2)

            plot_bins = np.arange(10000, 18000, 50)
            if grism.value == 'G102':
                plot_bins = np.arange(6000, 13000, 50)
            plot_spectrum = np.array([used_extraction_method(fits, ff, ff + 50, aperture_lower_extend.value,
                                                             aperture_upper_extend.value,
                                                             extraction_gauss_sigma.value)[0]
                                      for ff in plot_bins])

            plt.plot((plot_bins + 25) / 10000.0, plot_spectrum / 1000000.0, 'r-', lw=2)

            plt.ylabel(r'$\mathrm{e}^{-} \, (\times 10^6)$', fontsize=20)
            plt.xlabel(r'$\lambda \, [\mu \mathrm{m}]$', fontsize=20)

    for i in [heliocentric_julian_date_array, spectrum_direction_array, sky_background_level_array, x_star_array,
              x_shift_error_array, y_star_array, y_shift_error_array, scan_length_array, scan_length_error_array,
              white_dictionary] + bins_dictionaries:
        i.to_dictionary(light_curve)

    if plot:
        return input_data, [plt.figure(1), plt.figure(2)]

    else:
        return input_data, light_curve


def split_photometry(input_data, white_lower_wavelength=None, white_upper_wavelength=None, bins_file=None,
               aperture_lower_extend=None, aperture_upper_extend=None, extraction_method=None,
               extraction_gauss_sigma=None, plot=False):

    # load pipeline and calibration variables to be used

    white_lower_wavelength = variables.white_lower_wavelength.custom(white_lower_wavelength)
    white_upper_wavelength = variables.white_upper_wavelength.custom(white_upper_wavelength)
    bins_file = variables.bins_file.custom(bins_file)
    aperture_lower_extend = variables.aperture_lower_extend.custom(aperture_lower_extend)
    aperture_upper_extend = variables.aperture_upper_extend.custom(aperture_upper_extend)
    extraction_method = variables.extraction_method.custom(extraction_method)
    extraction_gauss_sigma = variables.extraction_gauss_sigma.custom(extraction_gauss_sigma)

    ra_target = variables.ra_target.custom()
    dec_target = variables.dec_target.custom()
    subarray_size = variables.sub_array_size.custom()
    grism = variables.grism.custom()
    exposure_time = variables.exposure_time.custom()
    bins_number = variables.bins_number.custom()
    heliocentric_julian_date_array = variables.heliocentric_julian_date_array.custom()
    spectrum_direction_array = variables.spectrum_direction_array.custom()
    sky_background_level_array = variables.sky_background_level_array.custom()
    x_star_array = variables.x_star_array.custom()
    x_shift_error_array = variables.x_shift_error_array.custom()
    y_star_array = variables.y_star_array.custom()
    y_shift_error_array = variables.y_shift_error_array.custom()
    scan_length_array = variables.scan_length_array.custom()
    scan_length_error_array = variables.scan_length_error_array.custom()
    white_ldc1 = variables.white_ldc1.custom()
    white_ldc2 = variables.white_ldc2.custom()
    white_ldc3 = variables.white_ldc3.custom()
    white_ldc4 = variables.white_ldc4.custom()

    lower_wavelength = variables.lower_wavelength.custom()
    upper_wavelength = variables.upper_wavelength.custom()
    flux_array = variables.flux_array.custom()
    error_array = variables.error_array.custom()
    ph_error_array = variables.ph_error_array.custom()

    # set bins

    white_dictionary, bins_dictionaries = \
        variables.set_binning(input_data, white_lower_wavelength.value, white_upper_wavelength.value,
                              white_ldc1.value, white_ldc2.value, white_ldc3.value, white_ldc4.value,
                              bins_file.value)

    # iterate over the splitted data sub-sets

    final_light_curve = {}

    for split_number, splitted_sub_set in enumerate(input_data.spectroscopic_images):

        if not plot:
            print('Splitting sample {0}:'.format(split_number + 1))

        light_curve = \
            photometry(input_data.copy_split(split_number),
                       white_lower_wavelength=white_lower_wavelength.value,
                       white_upper_wavelength=white_upper_wavelength.value,
                       aperture_lower_extend=aperture_lower_extend.value,
                       aperture_upper_extend=aperture_upper_extend.value,
                       bins_file=bins_file.value,
                       extraction_method=extraction_method.value,
                       extraction_gauss_sigma=extraction_gauss_sigma.value,
                       plot=False)[1]

        final_light_curve[variables.light_curve_split.keyword + str(split_number + 1)] = light_curve

        if plot:

            fits = splitted_sub_set[0]
            used_extraction_method = {'integral': get_flux_integral, 'gauss': get_flux_gauss}[extraction_method.value]

            total_plots = len(input_data.spectroscopic_images)
            plot_columns = 3
            plot_rows = int(total_plots / float(plot_columns) - 0.0000000001) + 1

            plt.figure(1, figsize=(3 * plot_columns, 3 * plot_rows))
            plt.subplot(plot_rows, plot_columns, split_number + 1)
            plt.title(r'{0}{1}{2}'.format('$\mathrm{split \, ', str(split_number + 1), '}$'), fontsize=20)
            if split_number + 1 != plot_columns * (plot_rows - 1) + 1:
                plt.tick_params(labelleft=False, labelbottom=False)

            plt.imshow(fits[1].data, origin='lower', aspect='auto')
            plt.xlim(0, len(fits[1].data))
            plt.ylim(0, len(fits[1].data))

            lower_wavelength.from_dictionary(white_dictionary)
            upper_wavelength.from_dictionary(white_dictionary)

            used_extraction_method(fits, lower_wavelength.value, upper_wavelength.value,
                                   aperture_lower_extend.value, aperture_upper_extend.value,
                                   extraction_gauss_sigma.value, plot=True)

            plt.figure(2)
            testx = np.arange(10000, 18000, 50)
            if grism.value == 'G102':
                testx = np.arange(6000, 13000, 50)
            testy = np.array([used_extraction_method(fits, ff, ff + 50, aperture_lower_extend.value,
                                                     aperture_upper_extend.value, extraction_gauss_sigma.value)[0]
                              for ff in testx])

            plt.plot((testx + 25) / 10000.0, testy / 1000000.0, '-', lw=2,
                     label=r'{0}{1}{2}'.format('$\mathrm{split \, ', str(split_number + 1), '}$'))

            if split_number + 1 == total_plots:

                plt.figure(1).text(0.005, 0.5, r'$\mathrm{row \, (pix)}$', fontsize=20,
                                   ha='center', va='center', rotation='vertical')
                plt.figure(1).text(0.5, 0.01, r'$\mathrm{column \, (pix)}$', fontsize=20, ha='center', va='center')

                plt.figure(2)
                plt.xlim(0.9, 2.0)
                plt.legend()
                plt.ylabel(r'$\mathrm{e}^{-} \, (\times 10^6)$', fontsize=20)
                plt.xlabel(r'$\lambda \, [\mu \mathrm{m}]$', fontsize=20)

        if split_number == 0:

            ra_target.from_dictionary(light_curve)
            ra_target.to_dictionary(final_light_curve)

            dec_target.from_dictionary(light_curve)
            dec_target.to_dictionary(final_light_curve)

            subarray_size.from_dictionary(light_curve)
            subarray_size.to_dictionary(final_light_curve)

            grism.from_dictionary(light_curve)
            grism.to_dictionary(final_light_curve)

            exposure_time.from_dictionary(light_curve)
            exposure_time.to_dictionary(final_light_curve)

            bins_number.from_dictionary(light_curve)
            bins_number.to_dictionary(final_light_curve)

            aperture_lower_extend.from_dictionary(light_curve)
            aperture_lower_extend.to_dictionary(final_light_curve)

            aperture_upper_extend.from_dictionary(light_curve)
            aperture_upper_extend.to_dictionary(final_light_curve)

            extraction_method.from_dictionary(light_curve)
            extraction_method.to_dictionary(final_light_curve)

            extraction_gauss_sigma.from_dictionary(light_curve)
            extraction_gauss_sigma.to_dictionary(final_light_curve)

            heliocentric_julian_date_array.from_dictionary(light_curve)
            heliocentric_julian_date_array.to_dictionary(final_light_curve)

            spectrum_direction_array.from_dictionary(light_curve)
            spectrum_direction_array.to_dictionary(final_light_curve)

            sky_background_level_array.from_dictionary(light_curve)
            sky_background_level_array.to_dictionary(final_light_curve)

            x_star_array.from_dictionary(light_curve)
            x_star_array.to_dictionary(final_light_curve)

            x_shift_error_array.from_dictionary(light_curve)
            x_shift_error_array.to_dictionary(final_light_curve)

            y_star_array.from_dictionary(light_curve)
            y_star_array.to_dictionary(final_light_curve)

            y_shift_error_array.from_dictionary(light_curve)
            y_shift_error_array.to_dictionary(final_light_curve)

            scan_length_array.from_dictionary(light_curve)
            scan_length_array.to_dictionary(final_light_curve)

            scan_length_error_array.from_dictionary(light_curve)
            scan_length_error_array.to_dictionary(final_light_curve)

            for i in [white_dictionary] + bins_dictionaries:
                i.from_dictionary(light_curve)
                i.to_dictionary(final_light_curve)

        else:

            for i in [white_dictionary] + bins_dictionaries:
                i.from_dictionary(light_curve)
                flux_array.from_dictionary(i)
                error_array.from_dictionary(i)
                ph_error_array.from_dictionary(i)
                current_flux = flux_array.value
                current_error = error_array.value
                current_ph_error = ph_error_array.value
                i.from_dictionary(final_light_curve)
                flux_array.from_dictionary(i)
                flux_array.to_dictionary(i, value=flux_array.value + current_flux)
                error_array.from_dictionary(i)
                error_array.to_dictionary(i, value=np.sqrt(error_array.value ** 2 + current_error ** 2))
                ph_error_array.from_dictionary(i)
                ph_error_array.to_dictionary(i, value=np.sqrt(ph_error_array.value ** 2 + current_ph_error ** 2))
                i.to_dictionary(final_light_curve)

    if plot:

        return input_data, [plt.figure(1), plt.figure(2)]

    else:

        return input_data, final_light_curve


def plot_photometry(dataset, lightcurve, directory):

    forward_colour = 'k'
    reverse_colour = 'r'

    hjd_time = lightcurve[variables.heliocentric_julian_date_array.keyword]
    flux = lightcurve[variables.white_dictionary.keyword][variables.flux_array.keyword]
    ssky = lightcurve[variables.sky_background_level_array.keyword]
    scan = lightcurve[variables.spectrum_direction_array.keyword]
    reverse = np.where(np.array(scan) < 0)
    forward = np.where(np.array(scan) > 0)

    if len(forward[0]) > 0:
        if dataset.splitted:
            original = list(dataset.spectroscopic_images)
            dataset.spectroscopic_images = list(np.array(dataset.spectroscopic_images)[:, forward[0]])
            test = []
            for i in dataset.spectroscopic_images:
                test.append([i[-1]])
            dataset.spectroscopic_images = test
            figures = split_photometry(dataset, plot=True)[1]
            dataset.spectroscopic_images = original
        else:
            if len(reverse[0]) > 0:
                save_the_reverse = dataset.spectroscopic_images[reverse[0][-1]]
            dataset.spectroscopic_images = [dataset.spectroscopic_images[forward[0][-1]]]
            figures = photometry(dataset, plot=True)[1]

        plc.save_figure(directory, figure=figures[0], name='forward_extraction_aperture')
        plc.save_figure(directory, figure=figures[1], name='forward_stellar_spectrum')

        plt.close('all')

    if len(reverse[0]) > 0:
        if dataset.splitted:
            dataset.spectroscopic_images = list(np.array(dataset.spectroscopic_images)[:, reverse[0]])
            test = []
            for i in dataset.spectroscopic_images:
                test.append([i[-1]])
            dataset.spectroscopic_images = test
            figures = split_photometry(dataset, plot=True)[1]
        else:
            dataset.spectroscopic_images = [save_the_reverse]
            figures = photometry(dataset, plot=True)[1]

        plc.save_figure(directory, figure=figures[0], name='reverse_extraction_aperture')
        plc.save_figure(directory, figure=figures[1], name='reverse_stellar_spectrum')

        plt.close('all')

    plt.subplot(2, 1, 1)
    plt.plot((np.array(hjd_time) - hjd_time[0])[forward], np.array(flux)[forward] / (10 ** 8),
             'o', c=forward_colour, mec=forward_colour, ms=3)
    plt.plot((np.array(hjd_time) - hjd_time[0])[reverse], np.array(flux)[reverse] / (10 ** 8),
             'o', c=reverse_colour, mec=reverse_colour, ms=3)
    plt.ylabel(r'$\mathrm{e}^{-} \, (\times 10^8)$', fontsize=15)
    plt.tick_params(labelbottom=False)
    plc.adjust_ticks()
    plt.subplot(2, 1, 2)
    plt.plot((np.array(hjd_time) - hjd_time[0])[forward], np.array(ssky)[forward],
             'o', c=forward_colour, mec=forward_colour, ms=3)
    plt.plot((np.array(hjd_time) - hjd_time[0])[reverse], np.array(ssky)[reverse],
             'o', c=reverse_colour, mec=reverse_colour, ms=3)
    plt.xlabel(r'$\Delta t \, \mathrm{(days)}$', fontsize=15)
    plt.ylabel(r'$\mathrm{sky} \, \mathrm{ratio}$', fontsize=15)
    plc.adjust_ticks()
    plt.subplots_adjust(hspace=0)
    plc.save_figure(directory, name='raw_light_curve')
