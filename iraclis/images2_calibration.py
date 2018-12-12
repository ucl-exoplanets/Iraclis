"""
images2_calibrationn.py

Includes all the functions that perform calibration processes.

All the main functions take as input either an HDUList object or a DataSet object, as defined in the basics.py file, and
return an updated version of the input. In all cases, the default values for the input parameters are the values in the
respective pipeline.variables object.

Main functions included:
calibration:                ...
split_recalibration:        ...

All the supporting functions take as input an HDUList object and return an updated version
(apart from the get_standard_flat function which does not affect the input). The parameters for the supporting functions
do not have default values, as their purpose is to be used only in this particular file.

Supporting functions included:
get_position_diagnostics:   ...
get_standard_flat:          ...
get_comparison_position:    ...
get_relative_position:      ...
get_scan_length:            ...
"""

from ._3objects import *


def get_position_diagnostics(fits):

    # set the reference pixels 0 so not to cause problems later

    first_read_index = functions.sci(fits)[-1]

    data = np.sum(fits[1].data, 0)
    model = tools.box(np.arange(len(data)), len(data) / 2, np.max(data), 45., 40.)
    dx = np.argmax(np.convolve(data, model)) - np.argmax(np.convolve(model, model))
    x_lim1, x_lim2 = int(len(data) / 2 + dx - 55), int(len(data) / 2 + dx + 55)

    # detect the approximate vertical position of the final spectrum

    data = fits[1].data
    y1 = 5 + np.median(np.argmax(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
    y2 = np.median(np.argmin(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
    if abs(y2 - y1) <= 1:
        y1 = 2 + np.median(np.argmax(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
        y2 = np.median(np.argmin(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
    final_y_lim1, final_y_lim2 = np.sort([int(round(float(y1))), int(round(float(y2)))])

    # detect the approximate vertical position of the spectrum in the first read

    data = fits[first_read_index].data
    y1 = 2.5 + np.median(np.argmax(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
    y2 = -2.5 + np.median(np.argmin(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
    if abs(y2 - y1) <= 1:
        y1 = 2 + np.median(np.argmax(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
        y2 = np.median(np.argmin(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
    initial_y_lim1, initial_y_lim2 = np.sort([int(round(float(y1))), int(round(float(y2)))])

    if abs(final_y_lim1 - final_y_lim2) > 1:
        final_scan = True
        if abs(initial_y_lim1 - initial_y_lim2) > 1:
            initial_scan = True
        else:
            initial_scan = False
        initial_y = 0.5 * (initial_y_lim1 + initial_y_lim2)
        direction = [1.0, -1.0][np.argmin(np.abs([final_y_lim1 - initial_y, final_y_lim2 - initial_y]))]

    else:
        final_scan = False
        initial_scan = False
        direction = 1.0

    return (final_y_lim1, final_y_lim2, x_lim1, x_lim2, final_scan, direction,
            initial_y_lim1, initial_y_lim2, initial_scan, direction)


def get_standard_flat(fits, sample, use_standard_flat):

    if use_standard_flat:

        spectrum_bottom = variables.spectrum_bottom.custom_from_fits(fits)
        spectrum_top = variables.spectrum_top.custom_from_fits(fits)
        spectrum_left = variables.spectrum_left.custom_from_fits(fits)
        spectrum_right = variables.spectrum_right.custom_from_fits(fits)

        flatfield = calibrations.flat_field_coefficient_1.match(fits)
        flatfield = np.where(flatfield == 0, 0.000001, flatfield)

        a = fits[1].data / flatfield
        detection_limit = 10

        ylim1 = spectrum_bottom.value - 5
        ylim2 = spectrum_top.value + 5
        xlim1 = spectrum_left.value
        xlim2 = spectrum_right.value

        xflag = np.abs(a - 0.25 * (np.roll(a, 1, 1) + np.roll(a, 2, 1) + np.roll(a, -1, 1) + np.roll(a, -2, 1)))
        medianx = np.median(xflag[ylim1:ylim2], 0)
        yflag = np.abs(a - 0.25 * (np.roll(a, 1, 0) + np.roll(a, 2, 0) + np.roll(a, -1, 0) + np.roll(a, -2, 0)))
        mediany = np.median(yflag[:, xlim1:xlim2], 1)
        medianx, mediany = np.meshgrid(medianx, mediany)
        madx = np.median(np.abs(xflag - medianx)[ylim1:ylim2], 0)
        mady = np.median(np.abs(yflag - mediany)[:, xlim1:xlim2], 1)
        madx, mady = np.meshgrid(madx, mady)
        limitx = medianx + detection_limit * madx
        limity = mediany + detection_limit * mady
        cr_test = np.where((xflag > limitx) & (yflag > limity), 1, 0)

        final_array = fits[sample].data / flatfield
        final_array = np.where(cr_test == 1, np.nan, final_array)

        while len(np.where(np.isnan(final_array))[0]) > 0:
            final_array = np.where(np.isnan(final_array),
                                   np.nanmean(
                                       [
                                           np.roll(final_array, 1, 0),
                                           np.roll(final_array, 2, 0),
                                           np.roll(final_array, -1, 0),
                                           np.roll(final_array, -2, 0),
                                           np.roll(final_array, 1, 1),
                                           np.roll(final_array, 2, 1),
                                           np.roll(final_array, -1, 1),
                                           np.roll(final_array, -2, 1),
                                       ], 0), final_array)

        return final_array

    else:

        return fits[sample].data


def get_absolute_x_star(fits, direct_image, target_x_offset):

    wfc3_aperture = variables.wfc3_aperture.custom()
    postarg1 = variables.postarg1.custom()
    grism = variables.grism.custom()

    direct_image_offsets = {'F098W': 0.150, 'F132N': 0.039, 'F140W': 0.083, 'F126N': 0.264, 'F153M': 0.146,
                            'F167N': 0.196, 'F139M': 0.110, 'F164N': 0.169, 'F127M': 0.131, 'F160W': 0.136,
                            'F128N': 0.026, 'F125W': 0.046, 'F130N': 0.033, 'F110W': -0.037, 'F105W': 0.015}

    direct_image_scales = {'IR-FIX': 0.135437, 'IR': 0.135601, 'G102-REF': 0.135603, 'G141-REF': 0.135603,
                           'IRSUB64': 0.135470, 'IRSUB128': 0.135470, 'IRSUB256': 0.135470, 'IRSUB512': 0.135470,
                           'IRSUB64-FIX': 0.135437, 'IRSUB128-FIX': 0.135437, 'IRSUB256-FIX': 0.135437,
                           'IRSUB512-FIX': 0.135437, 'IR-UVIS-CENTER': 0.135357, 'IR-UVIS': 0.135666,
                           'IR-UVIS-FIX': 0.135666, 'GRISM1024': 0.135603, 'GRISM512': 0.135504, 'GRISM256': 0.135508,
                           'GRISM128': 0.135404, 'GRISM64': 0.135404}

    direct_image_reference_pixel = {'IR-FIX': 512.0, 'IR': 562.0, 'G102-REF': 497.0, 'G141-REF': 497.0,
                                    'IRSUB64': 522.0, 'IRSUB128': 522.0, 'IRSUB256': 522.0, 'IRSUB512': 522.0,
                                    'IRSUB64-FIX': 512.0, 'IRSUB128-FIX': 512.0, 'IRSUB256-FIX': 512.0,
                                    'IRSUB512-FIX': 512.0, 'IR-UVIS-CENTER': 498.4, 'IR-UVIS': 490.9,
                                    'IR-UVIS-FIX': 490.9, 'GRISM1024': 497.0, 'GRISM512': 505.0, 'GRISM256': 410.0,
                                    'GRISM128': 496.0, 'GRISM64': 496.0}

    if fits[0].header['FILTER'] == 'G102':

        reference_pixel = {'IR-FIX': 512.0, 'IR': 497.0, 'G102-REF': 497.0, 'G141-REF': 497.0, 'IRSUB64': 522.0,
                           'IRSUB128': 522.0, 'IRSUB256': 522.0, 'IRSUB512': 522.0, 'IRSUB64-FIX': 512.0,
                           'IRSUB128-FIX': 512.0, 'IRSUB256-FIX': 512.0, 'IRSUB512-FIX': 512.0, 'IR-UVIS-CENTER': 498.4,
                           'IR-UVIS': 490.9, 'IR-UVIS-FIX': 490.9, 'GRISM1024': 497.0, 'GRISM512': 505.0,
                           'GRISM256': 410.0, 'GRISM128': 376.0, 'GRISM64': 376.0}

        scales = {'IR-FIX': 0.135437, 'IR': 0.135603, 'G102-REF': 0.135603, 'G141-REF': 0.135603, 'IRSUB64': 0.135470,
                  'IRSUB128': 0.135470, 'IRSUB256': 0.135470, 'IRSUB512': 0.135470, 'IRSUB64-FIX': 0.135437,
                  'IRSUB128-FIX': 0.135437, 'IRSUB256-FIX': 0.135437, 'IRSUB512-FIX': 0.135437,
                  'IR-UVIS-CENTER': 0.135357, 'IR-UVIS': 0.135666, 'IR-UVIS-FIX': 0.135666, 'GRISM1024': 0.135603,
                  'GRISM512': 0.135504, 'GRISM256': 0.135508, 'GRISM128': 0.135476, 'GRISM64': 0.135476}

    else:

        reference_pixel = {'IR-FIX': 512.0, 'IR': 497.0, 'G102-REF': 497.0, 'G141-REF': 497.0, 'IRSUB64': 522.0,
                           'IRSUB128': 522.0, 'IRSUB256': 522.0, 'IRSUB512': 522.0, 'IRSUB64-FIX': 512.0,
                           'IRSUB128-FIX': 512.0, 'IRSUB256-FIX': 512.0, 'IRSUB512-FIX': 512.0, 'IR-UVIS-CENTER': 498.4,
                           'IR-UVIS': 490.9, 'IR-UVIS-FIX': 490.9, 'GRISM1024': 497.0, 'GRISM512': 505.0,
                           'GRISM256': 410.0, 'GRISM128': 410.0, 'GRISM64': 410.0}

        scales = {'IR-FIX': 0.135437, 'IR': 0.135603, 'G102-REF': 0.135603, 'G141-REF': 0.135603, 'IRSUB64': 0.135470,
                  'IRSUB128': 0.135470, 'IRSUB256': 0.135470, 'IRSUB512': 0.135470, 'IRSUB64-FIX': 0.135437,
                  'IRSUB128-FIX': 0.135437, 'IRSUB256-FIX': 0.135437, 'IRSUB512-FIX': 0.135437,
                  'IR-UVIS-CENTER': 0.135357, 'IR-UVIS': 0.135666, 'IR-UVIS-FIX': 0.135666, 'GRISM1024': 0.135603,
                  'GRISM512': 0.135504, 'GRISM256': 0.135508, 'GRISM128': 0.135474, 'GRISM64': 0.135474}

    x0 = tools.fit_2d_gauss(direct_image[1].data)[0]

    grism.from_fits(direct_image)

    subarray_correction = 507 - len(direct_image[1].data[0]) / 2

    dx_off = direct_image_offsets[grism.value] - direct_image_offsets['F140W']

    wfc3_aperture.from_fits(direct_image)
    postarg1.from_fits(direct_image)

    x_ref_0 = (direct_image_reference_pixel[wfc3_aperture.value] + postarg1.value /
               direct_image_scales[wfc3_aperture.value])

    wfc3_aperture.from_fits(fits)
    postarg1.from_fits(fits)

    x_ref_1 = reference_pixel[wfc3_aperture.value] + postarg1.value / scales[wfc3_aperture.value]

    dx_ref = (x_ref_1 - x_ref_0)

    return x0 + subarray_correction + dx_off + dx_ref + target_x_offset


def get_absolute_y_star(fits, target_y_offset, use_standard_flat):

    spectrum_left = variables.spectrum_left.custom_from_fits(fits).value
    spectrum_right = variables.spectrum_right.custom_from_fits(fits).value
    first_spectrum_bottom = variables.first_spectrum_bottom.custom_from_fits(fits).value
    first_spectrum_top = variables.first_spectrum_top.custom_from_fits(fits).value
    spectrum_scan = variables.spectrum_scan.custom_from_fits(fits).value
    first_spectrum_scan = variables.first_spectrum_scan.custom_from_fits(fits).value
    first_spectrum_direction = variables.first_spectrum_direction.custom_from_fits(fits).value
    x_star = variables.x_star.custom_from_fits(fits).value

    trace_at0 = calibrations.trace_at0.match(fits)
    trace_at1 = calibrations.trace_at1.match(fits)
    trace_at2 = calibrations.trace_at2.match(fits)
    trace_at3 = calibrations.trace_at3.match(fits)
    trace_at4 = calibrations.trace_at4.match(fits)
    trace_at5 = calibrations.trace_at5.match(fits)
    trace_bt0 = calibrations.trace_bt0.match(fits)
    trace_bt1 = calibrations.trace_bt1.match(fits)
    trace_bt2 = calibrations.trace_bt2.match(fits)

    data = get_standard_flat(fits, functions.sci(fits)[-int(spectrum_scan)], use_standard_flat)
    data = np.swapaxes(data[max(5, first_spectrum_bottom - 20):min(first_spectrum_top + 20, len(fits[1].data) - 5),
                       spectrum_left:spectrum_right], 0, 1)
    rows = np.arange(max(5, first_spectrum_bottom - 20), min(first_spectrum_top + 20, len(fits[1].data) - 5))

    function_to_fit = [tools.fit_gauss, tools.fit_box][int(first_spectrum_scan)]
    avg = np.array([function_to_fit(rows, ff)[0] for ff in data])
    cols = np.arange(spectrum_left, spectrum_right)

    def trace(xarr, yy0):

        xx0 = x_star

        bt = trace_bt0 + trace_bt1 * xx0 + trace_bt2 * yy0
        at = (trace_at0 + trace_at1 * xx0 + trace_at2 * yy0 + trace_at3 * xx0 * xx0 +
              trace_at4 * xx0 * yy0 + trace_at5 * yy0 * yy0)
        subb = (507 - len(fits[1].data) / 2)

        return yy0 + bt + at * (xarr + subb - xx0) - subb

    y0 = curve_fit(trace, cols + 0.5, avg + 0.5, p0=[500])[0][0]

    if first_spectrum_scan:
        data = get_standard_flat(fits, functions.sci(fits)[-1], use_standard_flat)
        data = np.sum(data[5:-5, max(5, spectrum_left - 20):min(spectrum_right + 20, len(fits[1].data[0]) - 5)], 1)
        rows = np.arange(5, len(fits[1].data) - 5)
        center, center_err, fwhm, fwhm_err, popt = tools.fit_box(rows, data)
        correction = - first_spectrum_direction * fwhm / 2
    else:
        correction = 0

    return y0 + correction + target_y_offset


def calibration(input_data, direct_image=None, comparison_forward=None, comparison_index_forward=None,
                comparison_reverse=None, comparison_index_reverse=None, target_x_offset=None, target_y_offset=None,
                use_standard_flat=None, splitting=False):
    """
    Calibration.

    Calculates the horizontal position of the direct image of the star on the detector by comparing

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    direct_image : HDUList
        non-dispersed image of the target star that corresponds to the input fits file or to the comparison frame
        if the input_data object is a DataSet object and the direct_image is set to None, the direct image included
        in the DataSet object will be automatically used

    comparison_forward : HDUList
        comparison frame for the forward scans
        if the input_data object is an HDUList object and the comparison_forward is set to None,
        the input_data itself will be used as a comparison frame for the forward scans
        if the input_data object is a DataSet object and the comparison_forward is set to None,
        an HDUList from the DataSet object, indicated by the comparison_index_forward parameter,
        will be used as a comparison frame a comparison frame fro the forward scans

    comparison_index_forward : int
        index of the comparison frame for the forward scans in the case a DataSet object is given as input

    comparison_reverse : HDUList
        comparison frame for the reverse scans
        if the input_data object is an HDUList object and the comparison_reverse is set to None,
        the input_data itself will be used as a comparison frame fro the reverse scans
        if the input_data object is a DataSet object and the comparison_reverse is set to None,
        an HDUList from the DataSet object, indicated by the comparison_index_reverse parameter,
        will be used as a comparison frame a comparison frame fro the reverse scans

    comparison_index_reverse : int
        index of the comparison frame for the reverse scans in the case a DataSet object is given as input

    target_x_offset : float
        offset imposed by the user on the automatically calculated horizontal position of the direct image of the star
        on the detector

    target_y_offset : float
        offset imposed by the user on the automatically calculated vertical position of the direct image of the star
        on the detector at the beginning of the scan

    use_standard_flat : bool
        apply the standard non-wavelength-dependent flat-field as a correction step during the position calibration,
        this will not affect the data, which should be flat-fielded properly using the wavelength function

    plot : bool
        include the diagnostic plots in the output of the function

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, calibrated for position

    """

    input_data = DataSet(input_data)

    # load pipeline and calibration variables to be used

    target_x_offset = variables.target_x_offset.custom(target_x_offset)
    target_y_offset = variables.target_y_offset.custom(target_y_offset)
    comparison_index_forward = variables.comparison_index_forward.custom(comparison_index_forward)
    comparison_index_reverse = variables.comparison_index_reverse.custom(comparison_index_reverse)
    use_standard_flat = variables.use_standard_flat.custom(use_standard_flat)

    spectrum_bottom = variables.spectrum_bottom.custom()
    spectrum_top = variables.spectrum_top.custom()
    spectrum_left = variables.spectrum_left.custom()
    spectrum_right = variables.spectrum_right.custom()
    spectrum_scan = variables.spectrum_scan.custom()
    spectrum_direction = variables.spectrum_direction.custom()
    first_spectrum_bottom = variables.first_spectrum_bottom.custom()
    first_spectrum_top = variables.first_spectrum_top.custom()
    first_spectrum_scan = variables.first_spectrum_scan.custom()
    first_spectrum_direction = variables.first_spectrum_direction.custom()
    comparison_x_star = variables.comparison_x_star.custom()
    x_star = variables.x_star.custom()
    x_shift = variables.x_shift.custom()
    x_shift_error = variables.x_shift_error.custom()
    comparison_y_star = variables.comparison_y_star.custom()
    y_star = variables.y_star.custom()
    y_shift = variables.y_shift.custom()
    y_shift_error = variables.y_shift_error.custom()
    scan_length = variables.scan_length.custom()
    scan_length_error = variables.scan_length_error.custom()
    wdpt_constant_coefficient_1 = variables.wdpt_constant_coefficient_1.custom()
    wdpt_constant_coefficient_2 = variables.wdpt_constant_coefficient_2.custom()
    wdpt_constant_coefficient_3 = variables.wdpt_constant_coefficient_3.custom()
    wdpt_slope_coefficient_1 = variables.wdpt_slope_coefficient_1.custom()
    wdpt_slope_coefficient_2 = variables.wdpt_slope_coefficient_2.custom()
    wdpt_slope_coefficient_3 = variables.wdpt_slope_coefficient_3.custom()
    scan_frame = variables.scan_frame.custom()
    wavelength_frame = variables.wavelength_frame.custom()
    normalised_wavelength_frame = variables.normalised_wavelength_frame.custom()

    trace_at0 = calibrations.trace_at0.match(input_data)
    trace_at1 = calibrations.trace_at1.match(input_data)
    trace_at2 = calibrations.trace_at2.match(input_data)
    trace_at3 = calibrations.trace_at3.match(input_data)
    trace_at4 = calibrations.trace_at4.match(input_data)
    trace_at5 = calibrations.trace_at5.match(input_data)
    trace_bt0 = calibrations.trace_bt0.match(input_data)
    trace_bt1 = calibrations.trace_bt1.match(input_data)
    trace_bt2 = calibrations.trace_bt2.match(input_data)
    wsol_aw0 = calibrations.wsol_aw0.match(input_data)
    wsol_aw1 = calibrations.wsol_aw1.match(input_data)
    wsol_aw2 = calibrations.wsol_aw2.match(input_data)
    wsol_aw3 = calibrations.wsol_aw3.match(input_data)
    wsol_aw4 = calibrations.wsol_aw4.match(input_data)
    wsol_aw5 = calibrations.wsol_aw5.match(input_data)
    wsol_bw0 = calibrations.wsol_bw0.match(input_data)
    wsol_bw1 = calibrations.wsol_bw1.match(input_data)
    wsol_bw2 = calibrations.wsol_bw2.match(input_data)
    flat_field_min_wavelength = calibrations.flat_field_min_wavelength.match(input_data)
    flat_field_max_wavelength = calibrations.flat_field_max_wavelength.match(input_data)

    if not splitting:

        # set up the direct image and the comparison scans

        if isinstance(input_data, DataSet):
            if direct_image is None:
                direct_image = functions.fits_like(input_data.direct_image)
            else:
                direct_image = functions.fits_like(direct_image)
            if comparison_forward is None:
                comparison_forward = functions.fits_like(input_data.spectroscopic_images[comparison_index_forward.value])
            else:
                comparison_forward = functions.fits_like(comparison_forward)
            if comparison_reverse is None:
                comparison_reverse = functions.fits_like(input_data.spectroscopic_images[comparison_index_reverse.value])
            else:
                comparison_reverse = functions.fits_like(comparison_reverse)
        else:
            if direct_image is None:
                raise IraclisInputError('Direct image not given.')
            else:
                direct_image = functions.fits_like(direct_image)
            if comparison_forward is None:
                comparison_forward = functions.fits_like(input_data)
            else:
                comparison_forward = functions.fits_like(comparison_forward)
            if comparison_reverse is None:
                comparison_reverse = functions.fits_like(input_data)
            else:
                comparison_reverse = functions.fits_like(comparison_reverse)

        comparisons = {'forward': {'x_star': None, 'y_star': None, 'spectrum_right': None, 'spectrum_bottom': None,
                                   'interp_x': None, 'interp_y': None},
                       'reverse': {'x_star': None, 'y_star': None, 'spectrum_right': None, 'spectrum_bottom': None,
                                   'interp_x': None, 'interp_y': None}}

        for name, comparison in [['forward', comparison_forward], ['reverse', comparison_reverse]]:

            for i in functions.sci(comparison):
                comparison[i].data[:5, :] = 0
                comparison[i].data[-5:, :] = 0
                comparison[i].data[:, :5] = 0
                comparison[i].data[:, -5:] = 0

            diagnostics = get_position_diagnostics(comparison)
            spectrum_bottom.to_fits(comparison, value=diagnostics[0])
            spectrum_top.to_fits(comparison, value=diagnostics[1])
            spectrum_left.to_fits(comparison, value=diagnostics[2])
            spectrum_right.to_fits(comparison, value=diagnostics[3])
            spectrum_scan.to_fits(comparison, value=diagnostics[4])
            spectrum_direction.to_fits(comparison, value=diagnostics[5])
            first_spectrum_bottom.to_fits(comparison, value=diagnostics[6])
            first_spectrum_top.to_fits(comparison, value=diagnostics[7])
            first_spectrum_scan.to_fits(comparison, value=diagnostics[8])
            first_spectrum_direction.to_fits(comparison, value=diagnostics[9])

            x_star.to_fits(comparison, value=get_absolute_x_star(comparison, direct_image, target_x_offset.value))
            y_star.to_fits(comparison, value=get_absolute_y_star(comparison, target_y_offset.value,
                                                                 use_standard_flat.value))

            comparisons[name]['x_star'] = x_star.value
            comparisons[name]['y_star'] = y_star.value
            comparisons[name]['spectrum_right'] = spectrum_right.value
            comparisons[name]['spectrum_bottom'] = first_spectrum_bottom.value

            data = get_standard_flat(comparison, functions.sci(comparison)[0], use_standard_flat.value)
            if spectrum_scan.value and abs(spectrum_top.value - spectrum_bottom.value) > 15:
                data = np.sum(data[spectrum_bottom.value + 5:spectrum_top.value - 5], 0)
            else:
                data = np.sum(data[spectrum_bottom.value - 5:spectrum_top.value + 5], 0)

            cols = np.arange(len(data))
            for expand in range(5, 200):
                cols = np.append(cols, len(comparison[1].data) + expand)
                data = np.append(data, 0)
                cols = np.append(cols, - expand)
                data = np.append(data, 0)

            comparisons[name]['interp_x'] = interp1d(cols, data, kind='cubic')

            data = get_standard_flat(comparison,
                                     functions.sci(comparison)[-int(spectrum_scan.value)], use_standard_flat.value)
            data = np.sum(
                data[:, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                             len(comparison[1].data) - 7)], 1)

            rows = np.arange(len(data))
            for expand in range(5, 200):
                rows = np.append(rows, len(comparison[1].data) + expand)
                data = np.append(data, 0)
                rows = np.append(rows, - expand)
                data = np.append(data, 0)

            comparisons[name]['interp_y'] = interp1d(rows, data, kind='cubic')

        # initiate counter

        counter = PipelineCounter('Calibration', len(input_data.spectroscopic_images))

        # iterate over the list of HDUList objects included in the input data

        for fits in input_data.spectroscopic_images:

            # run the position diagnostics

            for i in functions.sci(fits):
                fits[i].data[:5, :] = 0
                fits[i].data[-5:, :] = 0
                fits[i].data[:, :5] = 0
                fits[i].data[:, -5:] = 0

            diagnostics = get_position_diagnostics(fits)
            spectrum_bottom.to_fits(fits, value=diagnostics[0])
            spectrum_top.to_fits(fits, value=diagnostics[1])
            spectrum_left.to_fits(fits, value=diagnostics[2])
            spectrum_right.to_fits(fits, value=diagnostics[3])
            spectrum_scan.to_fits(fits, value=diagnostics[4])
            spectrum_direction.to_fits(fits, value=diagnostics[5])
            first_spectrum_bottom.to_fits(fits, value=diagnostics[6])
            first_spectrum_top.to_fits(fits, value=diagnostics[7])
            first_spectrum_scan.to_fits(fits, value=diagnostics[8])
            first_spectrum_direction.to_fits(fits, value=diagnostics[9])

            # choose the comparison scan based on the scan direction

            if spectrum_direction.value > 0:
                comparison = comparison_forward
                comparison_name = 'forward'
            else:
                comparison = comparison_reverse
                comparison_name = 'reverse'

            # calculate the horizontal shift

            testx = np.arange(len(fits[1].data))
            testy = get_standard_flat(fits, functions.sci(fits)[0], use_standard_flat.value)

            spectrum_bottom.from_fits(comparison)
            spectrum_top.from_fits(comparison)
            spectrum_scan.from_fits(comparison)

            if spectrum_scan.value and abs(spectrum_top.value - spectrum_bottom.value) > 15:
                testy = np.sum(testy[spectrum_bottom.value + 5:spectrum_top.value - 5], 0)
            else:
                testy = np.sum(testy[spectrum_bottom.value - 5:spectrum_top.value + 5], 0)

            spectrum_bottom.from_fits(fits)
            spectrum_top.from_fits(fits)
            spectrum_scan.from_fits(fits)

            # plt.plot(testy)

            interp_x = comparisons[comparison_name]['interp_x']

            def fit_xshift(xx, ddx, ddf):
                return ddf * interp_x(xx + ddx)

            margin1 = max(spectrum_left.value - 40, 6)
            margin2 = min(spectrum_right.value + 40, len(fits[1].data) - 6)
            best_fit, covariance = curve_fit(fit_xshift, testx[margin1:margin2], testy[margin1:margin2],
                                             p0=[-spectrum_right.value + comparisons[comparison_name]['spectrum_right'],
                                             np.median(testy[margin1:margin2] / interp_x(testx[margin1:margin2]))])

            comparison_x_star.to_fits(fits, value=comparisons[comparison_name]['x_star'])
            x_star.to_fits(fits, value=comparisons[comparison_name]['x_star'] - best_fit[0])
            x_shift.to_fits(fits, value=-best_fit[0])
            x_shift_error.to_fits(fits, value=[np.sqrt(covariance[0][0]), 0][int(np.isinf(np.sqrt(covariance[0][0])))])

            # calculate the vertical shift

            testx = np.arange(len(fits[1].data))
            testy = get_standard_flat(fits, functions.sci(fits)[-int(spectrum_scan.value)], use_standard_flat.value)

            testy = np.sum(
                testy[:, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                              len(comparison[1].data) - 7)], 1)

            interp_y = comparisons[comparison_name]['interp_y']

            def fit_yshift(yyy, ddy, ddf):
                return ddf * interp_y(yyy + ddy)

            margin1 = max(first_spectrum_bottom.value - 20, 6)
            margin2 = min(first_spectrum_top.value + 20, len(fits[1].data) - 6)
            best_fit, covariance = curve_fit(fit_yshift, testx[margin1:margin2], testy[margin1:margin2],
                                             p0=[-first_spectrum_bottom.value +
                                                 comparisons[comparison_name]['spectrum_bottom'],
                                             np.median(testy[margin1:margin2] / interp_y(testx[margin1:margin2]))])

            comparison_y_star.to_fits(fits, value=comparisons[comparison_name]['y_star'])
            y_star.to_fits(fits, value=comparisons[comparison_name]['y_star'] - best_fit[0])
            y_shift.to_fits(fits, value=-best_fit[0])
            y_shift_error.to_fits(fits, value=[np.sqrt(covariance[0][0]), 0][int(np.isinf(np.sqrt(covariance[0][0])))])

            # calculate the scan length

            if target_x_offset.value != 0 or target_y_offset.value != 0:
                data = get_standard_flat(fits, 1, use_standard_flat.value)
                data = np.sum(data[spectrum_bottom.value - 10:spectrum_top.value + 10,
                              max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                   len(fits[1].data) - 7)], 1)
                rows = np.arange(spectrum_bottom.value - 10, spectrum_top.value + 10)
                for i in range(100):
                    data = np.append(data, 0)
                    rows = np.append(rows, spectrum_top.value + 20 + i)
                    data = np.append(data, 0)
                    rows = np.append(rows, spectrum_bottom.value - 20 - i)
            else:
                data = get_standard_flat(fits, 1, use_standard_flat.value)
                data = np.sum(
                    data[5:-5, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                    len(fits[1].data) - 7)], 1)
                rows = np.arange(5, len(fits[1].data) - 5)

            if spectrum_scan.value:
                best_fit = tools.fit_box(rows, data)[2:4]
            else:
                best_fit = [0, 0]

            scan_length.to_fits(fits, value=best_fit[0])
            scan_length_error.to_fits(fits, value=best_fit[1])

            # wavelength calibration

            y0 = np.arange(507 - len(fits[1].data) / 2, 507 + len(fits[1].data) / 2, 10)
            wl = np.arange(10000, 19000, 10)
            y0, wl = np.meshgrid(y0, wl)
            x0 = np.ones_like(y0) * x_star.value

            bt = trace_bt0 + trace_bt1 * x0 + trace_bt2 * y0
            at = (trace_at0 + trace_at1 * x0 + trace_at2 * y0 + trace_at3 * x0 * x0 +
                  trace_at4 * x0 * y0 + trace_at5 * y0 * y0)
            bw = wsol_bw0 + wsol_bw1 * x0 + wsol_bw2 * y0
            aw = (wsol_aw0 + wsol_aw1 * x0 + wsol_aw2 * y0 + wsol_aw3 * x0 * x0 +
                  wsol_aw4 * x0 * y0 + wsol_aw5 * y0 * y0)

            xl = x0 - (bt * at) / (at ** 2 + 1.0) + ((wl - bw) / aw) * np.cos(np.arctan(at))
            yl = at * (xl - x0) + bt + y0

            def level(xy, a, b, cc, d, e, f):
                x_lev, y_lev = xy
                return a * x_lev + b * y_lev + cc * x_lev * x_lev + d * x_lev * y_lev + e * y_lev * y_lev + f

            def wdpt(xlam, c1, c2, c3, s1, s2, s3):
                x_wdpt, l_wdpt = xlam
                return (c1 / (c2 + l_wdpt) + c3) + (s1 / (s2 + l_wdpt) + s3) * x_wdpt

            cor = (507 - len(fits[1].data) / 2)

            best_fit, covariance = curve_fit(level, (xl.flatten() - cor, yl.flatten() - cor), wl.flatten(),
                                             p0=[40, 1, 1, 1, 1, -1000])

            wlcoeff1, wlcoeff2, wlcoeff3, wlcoeff4, wlcoeff5, wlcoeff6 = best_fit

            best_fit, covariance = curve_fit(level, (xl.flatten() - cor, yl.flatten() - cor), y0.flatten(),
                                             p0=[40, 1, 1, 1, 1, -1000])

            y0coeff1, y0coeff2, y0coeff3, y0coeff4, y0coeff5, y0coeff6 = best_fit

            try:
                tests = []
                for ap in [-10 ** ff for ff in range(6, 9)] + [10 ** ff for ff in range(6, 9)]:
                    for bp in [1, 10000, -10000]:
                        for cp in [10000, -10000]:
                            for dp in [-10 ** 6, 10 ** 6]:
                                for ep in [1, -10 ** 4, 10 ** 4]:
                                    test = np.sum((wdpt((xl.flatten() - cor, wl.flatten()), ap, bp, cp, dp, ep, 1) -
                                                   (yl.flatten() - cor))**2)
                                    tests.append([test, [ap, bp, cp, dp, ep, 1]])

                tests.sort()
                best_fit, covariance = curve_fit(wdpt, (xl.flatten() - cor, wl.flatten()), yl.flatten() - cor,
                                                 p0=tests[0][1])
            except RuntimeError:
                if len(fits[1].data) < 200:
                    best_fit, covariance = curve_fit(wdpt, (xl.flatten() - cor, wl.flatten()), yl.flatten() - cor,
                                                     p0=[-10 ** 8, 1, -10000, -10 ** 6, 1, 1], maxfev=20000)
                elif len(fits[1].data) < 400:
                    best_fit, covariance = curve_fit(wdpt, (xl.flatten() - cor, wl.flatten()), yl.flatten() - cor,
                                                     p0=[10 ** 6, 1, -10000, -10 ** 6, 1, 1], maxfev=20000)
                else:
                    if x_star.value > 450:
                        best_fit, covariance = curve_fit(wdpt, (xl.flatten() - cor, wl.flatten()), yl.flatten() - cor,
                                                         p0=[10 ** 8, 1, -10000, -10 ** 6, 1, 1], maxfev=20000)
                    else:
                        best_fit, covariance = curve_fit(wdpt, (xl.flatten() - cor, wl.flatten()), yl.flatten() - cor,
                                                         p0=[-10 ** 8, 1, -10000, -10 ** 6, 1, 1], maxfev=20000)

            wdpt_constant_coefficient_1.to_fits(fits, value=best_fit[0])
            wdpt_constant_coefficient_2.to_fits(fits, value=best_fit[1])
            wdpt_constant_coefficient_3.to_fits(fits, value=best_fit[2])
            wdpt_slope_coefficient_1.to_fits(fits, value=best_fit[3])
            wdpt_slope_coefficient_2.to_fits(fits, value=best_fit[4])
            wdpt_slope_coefficient_3.to_fits(fits, value=best_fit[5])

            x, y = np.meshgrid(np.arange(len(fits[1].data)) + 0.5, np.arange(len(fits[1].data)) + 0.5)
            yy = (y0coeff1 * x + y0coeff2 * y + y0coeff3 * x * x + y0coeff4 * x * y + y0coeff5 * y * y + y0coeff6)
            scan_frame.to_fits(fits, value=yy - y_star.value)
            ww = (wlcoeff1 * x + wlcoeff2 * y + wlcoeff3 * x * x + wlcoeff4 * x * y + wlcoeff5 * y * y + wlcoeff6)
            wavelength_frame.to_fits(fits, value=ww)

            wmin = flat_field_min_wavelength
            wmax = flat_field_max_wavelength
            ww = np.where(scan_frame.value < - 100, wmin, ww)
            ww = np.where(scan_frame.value > scan_length.value + 100, wmin, ww)
            ww = np.where(ww < wmin, wmin, ww)
            ww = np.where(ww > wmax, wmin, ww)
            ww = (ww - wmin) / (wmax - wmin)
            normalised_wavelength_frame.to_fits(fits, value=ww)

            # update counter

            counter.update()

    else:

        comparison_forward = functions.fits_like(input_data.spectroscopic_images[comparison_index_forward.value])
        comparison_reverse = functions.fits_like(input_data.spectroscopic_images[comparison_index_reverse.value])

        # calibrate comparisons for vertical position

        comparisons = {'forward': {'y_star': None, 'spectrum_bottom': None, 'interp_y': None, 'scan': None},
                       'reverse': {'y_star': None, 'spectrum_bottom': None, 'interp_y': None, 'scan': None}}

        for name, comparison in [['forward', comparison_forward], ['reverse', comparison_reverse]]:

            for i in functions.sci(comparison):
                comparison[i].data[:5, :] = 0
                comparison[i].data[-5:, :] = 0
                comparison[i].data[:, :5] = 0
                comparison[i].data[:, -5:] = 0

            diagnostics = get_position_diagnostics(comparison)
            spectrum_bottom.to_fits(comparison, value=diagnostics[0])
            spectrum_top.to_fits(comparison, value=diagnostics[1])
            spectrum_scan.to_fits(comparison, value=diagnostics[4])
            first_spectrum_bottom.to_fits(comparison, value=diagnostics[6])
            first_spectrum_top.to_fits(comparison, value=diagnostics[7])
            first_spectrum_scan.to_fits(comparison, value=diagnostics[8])

            spectrum_left.from_fits(comparison)
            spectrum_right.from_fits(comparison)
            x_star.from_fits(comparison)
            y_star.to_fits(comparison, value=get_absolute_y_star(comparison, target_y_offset.value,
                                                                 use_standard_flat.value))

            comparisons[name]['scan'] = spectrum_scan.value
            comparisons[name]['y_star'] = y_star.value
            comparisons[name]['spectrum_bottom'] = spectrum_bottom.value

            data = get_standard_flat(comparison,
                                     functions.sci(comparison)[-int(spectrum_scan.value)], use_standard_flat.value)
            data = np.sum(data[:, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                       len(comparison[1].data) - 7)], 1)

            rows = np.arange(len(data))
            for expand in range(5, 200):
                rows = np.append(rows, len(comparison[1].data) + expand)
                data = np.append(data, 0)
                rows = np.append(rows, - expand)
                data = np.append(data, 0)

            comparisons[name]['interp_y'] = interp1d(rows, data, kind='cubic')

        # initiate counter

        counter = PipelineCounter('Calibration', len(input_data.spectroscopic_images))

        # iterate over the list of HDUList objects included in the input data

        for fits in input_data.spectroscopic_images:

            for i in functions.sci(fits):
                fits[i].data[:5, :] = 0
                fits[i].data[-5:, :] = 0
                fits[i].data[:, :5] = 0
                fits[i].data[:, -5:] = 0

            # run the position diagnostics

            diagnostics = get_position_diagnostics(fits)
            spectrum_bottom.to_fits(fits, value=diagnostics[0])
            spectrum_top.to_fits(fits, value=diagnostics[1])
            first_spectrum_bottom.to_fits(fits, value=diagnostics[6])
            first_spectrum_top.to_fits(fits, value=diagnostics[7])

            spectrum_direction.from_fits(fits)
            spectrum_left.from_fits(fits)
            spectrum_right.from_fits(fits)

            y_star.from_fits(fits)
            original_star_y_position = y_star.value
            scan_frame.from_fits(fits)
            original_scan_frame = scan_frame.value

            # choose the comparison scan based on the scan direction

            if spectrum_direction.value > 0:
                comparison = comparison_forward
                comparison_name = 'forward'
            else:
                comparison = comparison_reverse
                comparison_name = 'reverse'

            # calculate the vertical shift

            spectrum_scan.to_fits(fits, value=comparisons[comparison_name]['scan'])
            first_spectrum_scan.to_fits(fits, value=comparisons[comparison_name]['scan'])

            testx = np.arange(len(fits[1].data))
            testy = get_standard_flat(fits, functions.sci(fits)[-int(spectrum_scan.value)], use_standard_flat.value)
            testy = np.sum(testy[:, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                         len(comparison[1].data) - 7)], 1)

            interp_y = comparisons[comparison_name]['interp_y']

            def fit_yshift(y, ddy, ddf):
                return ddf * interp_y(y + ddy)

            margin1 = max(first_spectrum_bottom.value - 20, 6)
            margin2 = min(first_spectrum_top.value + 20, len(fits[1].data) - 6)
            best_fit, covariance = curve_fit(fit_yshift, testx[margin1:margin2], testy[margin1:margin2],
                                             p0=[-spectrum_bottom.value +
                                                 comparisons[comparison_name]['spectrum_bottom'],
                                             np.median(testy[margin1:margin2] / interp_y(testx[margin1:margin2]))])

            comparison_y_star.to_fits(fits, value=comparisons[comparison_name]['y_star'])
            y_star.to_fits(fits, value=comparisons[comparison_name]['y_star'] - best_fit[0])
            y_shift.to_fits(fits, value=-best_fit[0])
            y_shift_error.to_fits(fits, value=[np.sqrt(covariance[0][0]), 0][int(np.isinf(np.sqrt(covariance[0][0])))])

            # calculate the scan length

            if target_y_offset.value != 0 or target_y_offset.value != 0:
                data = get_standard_flat(fits, 1, use_standard_flat.value)
                data = np.sum(data[spectrum_bottom.value - 10:spectrum_top.value + 10,
                              max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                   len(fits[1].data) - 7)], 1)
                rows = np.arange(spectrum_bottom.value - 10, spectrum_top.value + 10)
                for i in range(100):
                    data = np.append(data, 0)
                    rows = np.append(rows, spectrum_top.value + 20 + i)
                    data = np.append(data, 0)
                    rows = np.append(rows, spectrum_bottom.value - 20 - i)
            else:
                data = get_standard_flat(fits, 1, use_standard_flat.value)
                data = np.sum(
                    data[5:-5, max(7, spectrum_left.value - 20):min(spectrum_right.value + 20,
                                                                    len(fits[1].data) - 7)], 1)
                rows = np.arange(5, len(fits[1].data) - 5)

            if spectrum_scan.value:
                best_fit = tools.fit_box(rows, data)[2:4]
            else:
                best_fit = [0, 0]

            scan_length.to_fits(fits, value=best_fit[0])
            scan_length_error.to_fits(fits, value=best_fit[1])

            scan_frame.to_fits(fits, value=original_scan_frame - (y_star.value - original_star_y_position))

            wavelength_frame.from_fits(fits)
            ww = wavelength_frame.value
            wmin = flat_field_min_wavelength
            wmax = flat_field_max_wavelength
            ww = np.where(scan_frame.value < - 100, wmin, ww)
            ww = np.where(scan_frame.value > scan_length.value + 100, wmin, ww)
            ww = np.where(ww < wmin, wmin, ww)
            ww = np.where(ww > wmax, wmin, ww)
            ww = (ww - wmin) / (wmax - wmin)
            normalised_wavelength_frame.to_fits(fits, value=ww)

            # update counter

            counter.update()

    # return the updated version of the input data

    return input_data
