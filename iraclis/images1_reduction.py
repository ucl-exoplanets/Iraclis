"""
File: images1_reduction.py
Author: Angelos Tsiaras
E-mail: angelos.tsiaras.14@ucl.ac.uk
Contents: functions that perform reduction processes
    timing:     Heliocentric julian date calculation.
    bias:       Bias drifts and zero-read corrections.
    linearity:
    dark:
    gain:       Gain variations correction and DN to electrons conversion.
    sky:
    flat:
    bpcr:

"""

from basics import *


def timing(input_data):
    """
    Heliocentric julian date calculation.

    Calculates the heliocentric julian date of the mid-exposure (hjd):

    .. math:: hjd = jd - (d / c) (\sin(d) \sin(ds_{jd}) + \cos(d) \cos(ds_{jd}) \cos(a - as_{jd})

    with:

    .. math:: jd = 0.5 (exp_s + exp_e) + 2400000.5

    where:

    jd is the julian date of the mid-exposure,
    d is the distance between the earth and the sun, 149597870700.0 km,
    c is the speed of light, 25902068371200.0 km/day,
    exp_s and exp_e are the exposure start and end times in modified julian date (included in the data file),
    a and d are the target RA and DEC (included in the data file) and
    as_jd and ds_jd are the sun RA and DEC at the julian date of the mid-exposure.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, including the heliocentric julian date of the mid-exposure

    """

    # load pipeline and calibration variables to be used

    exposure_start = pipeline_variables.exposure_start.custom()
    exposure_end = pipeline_variables.exposure_end.custom()
    ra_target = pipeline_variables.ra_target.custom()
    dec_target = pipeline_variables.dec_target.custom()
    heliocentric_julian_date = pipeline_variables.heliocentric_julian_date.custom()

    # initiate counter

    counter = PipelineCounter('Timing', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # get the exposure start and end times, as well as the target RA and DEC from the fits file

        exposure_start.from_fits(fits)
        exposure_end.from_fits(fits)
        ra_target.from_fits(fits)
        dec_target.from_fits(fits)

        # calculate the julian date of the mid-exposure

        julian_date = 0.5 * (exposure_start.value + exposure_end.value) + 2400000.5

        # calculate the sun RA and DEC at the julian date of the mid-exposure

        sun = ephem.Sun()
        sun.compute(ephem.date(julian_date - 2415020.0))  # the ephem package uses dublin julian date (-2415020 days)
        ra_sun_jd, dec_sun_jd = float(sun.ra), float(sun.dec)

        # calculate (and save in the data file) the heliocentric julian date of the mid-exposure, in the data files the
        # target RA and DEC are provided in degrees, hence we need to multiply them by np.pi / 180

        heliocentric_julian_date.to_fits(fits, value=(julian_date - ((149597870700.0 / 25902068371200.0) *
                                                      (np.sin(dec_target.value * np.pi / 180) * np.sin(dec_sun_jd) +
                                                       np.cos(dec_target.value * np.pi / 180) * np.cos(dec_sun_jd) *
                                                       np.cos(ra_target.value * np.pi / 180 - ra_sun_jd)))))

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def bias(input_data):
    """
    Bias drifts and zero-read corrections.

    1) Corrects each NDR for the zero-read.

    As the WFC3/IR detector lacks a shutter, the pixels are exposed even before the exposure starts. Hence the ZR is a
    record of the light collected before the beginning of the exposure. Since this light is not part of the exposure,
    the ZR is subtracted from each consecutive NDR:

    .. math:: NDR \rightarrow NDR - ZR

    At this stage the first NDR is renamed from SCI, not to be considered as an NDR for the rest of the analysis.

    2) Corrects each NDR for the bias drifts.

    The reference pixels define a group of pixels that are not sensitive to incoming light and are located at the first
    five and last five rows and columns of the detector. If the bias level is varying between the zero-read and an NDR,
    the variation is imprinted in the reference pixels. For this reason, the mean value of the reference pixels
    (reference level, NDR_rl) is subtracted from each NDR. The first five and the last five rows, as well as the
    first and the last columns, are excluded from the calculation because they have been found to be unstable.

    .. math:: NDR \rightarrow NDR - NDR_{rl}

    3) Calculates the error-array of each NDR.

    This initial calculation of the error-array of an NDR (NDR_e), combines the detector read noise (RN) and the
    Poison noise. The detector read noise is included in the WFC3 ccd calibration file in electron units, while the
    Poison noise is the square root of the number of electrons recorded in an NDR. At this stage the NDR is in DN units,
    so the number of electrons is given by multiplying the NDR by the detector gain (G), which is also included in the
    WFC3 ccd calibration file. The final uncertainty is converted back to DN units by dividing by the detector gain.

    .. math:: NDR_e (e- units) = \sqrt( RN ^ 2 +  NDR * G )
    .. math:: NDR_e (DN units) = \sqrt( RN ^ 2 +  NDR * G ) / G

    The error-array of each NDR is stored in the data file just after the NDR.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for bias drifts and zero-read

    """

    # load pipeline and calibration variables to be used

    zero_read_frame = pipeline_variables.zero_read.custom()
    zero_read_error_frame = pipeline_variables.zero_read_error.custom()
    reference_pixels_level = pipeline_variables.reference_pixels_level.custom()
    
    ccd_gain = calibration_variables.ccd_gain.match(input_data)
    ccd_read_noise = calibration_variables.ccd_read_noise.match(input_data)

    # initiate counter

    counter = PipelineCounter('Bias', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # check if the zero-read has been identified

        try:

            # get the zero-read and the zero-read error-array from the fits file

            zero_read_frame.from_fits(fits)
            zero_read_error_frame.from_fits(fits)

        except KeyError:

            # locate and rename the zero-read and the zero-read error-array
            # (first non-destructive read - i.e. last in the list of scientific frames)

            fits[sci(fits)[-1]].name = zero_read_frame.keyword
            fits[err(fits)[-1]].name = zero_read_error_frame.keyword

            # get the zero-read and the zero-read error-array from the fits file

            zero_read_frame.from_fits(fits)
            zero_read_error_frame.from_fits(fits)

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR from the fits file

            science = np.array(fits[i].data, dtype=float)

            # correct the NDR for the the zero-read

            science -= zero_read_frame.value

            # calculate (and save in the data file) the reference level

            reference_pixels_level.to_fits(
                fits,
                value=np.mean(np.sort(np.insert(science[5:-5, -5:-1].flatten(), 0, science[5:-5, 1:5].flatten()))),
                position=i)
            
            # correct the NDR for the the bias drifts

            science -= reference_pixels_level.value

            # calculate the error-array of the NDR

            error = (np.sqrt(ccd_read_noise * ccd_read_noise + ccd_gain * np.abs(science)) / ccd_gain)

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def linearity(input_data):
    """
    Non-linearity correction.

    1) Calculates the zero-read flux.

    The first NDR (zero-read, ZR) is a combination of the default bias level of the detector (super-zero-read, SZR) and
    the flux recorded before the beginning of the exposure (zero-read flux, ZRF). The super-zero-read is included in the
    WFC3 linearity calibration file. Given the above, the zero-read flux is calculated by subtracting from the zero-read
    the super-zero-read:

    .. math:: ZRF = ZR - SZR

    The ZRF is stored in the data file as the last element of the HDUList object.





    2) Calculates the error-array of the zero-read flux.



    .. math:: ZRF_e = \sqrt{RN ^ 2 + (SZR_e * G ) ^ 2 + ZRF * G} / G







    3) Corrects each read for the non-linearity behaviour of the IR detector, by applying on it the function:

    .. math:: F_{final} = F_c(f_r + f_z) - F_c(f_z)

    where fr is the flux in the read, fz is the zero-read flux, and Fc is the non-linearity correction function:

    .. math:: F_c(f) = (1 + c1 + c2 * f + c3 * f ^ 2 + c4 * f ^ 3) * f

    where c1, c2, c3, c4 are the non-linearity coefficient arrays included in the linearity calibration file.



    4)

    Also propagates the uncertainties in the error arrays.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for non-linearity

    """

    # load pipeline and calibration variables to be used

    zero_read = pipeline_variables.zero_read.custom()
    zero_read_error = pipeline_variables.zero_read_error.custom()
    zero_read_flux = pipeline_variables.zero_read_flux.custom()
    zero_read_flux_error = pipeline_variables.zero_read_flux_error.custom()

    super_zero_read = calibration_variables.super_zero_read.match(input_data)
    super_zero_read_error = calibration_variables.super_zero_read_error.match(input_data)
    ccd_gain = calibration_variables.ccd_gain.match(input_data)
    ccd_read_noise = calibration_variables.ccd_read_noise.match(input_data)
    linearity_coefficient_1 = calibration_variables.linearity_coefficient_1.match(input_data)
    linearity_coefficient_2 = calibration_variables.linearity_coefficient_2.match(input_data)
    linearity_coefficient_3 = calibration_variables.linearity_coefficient_3.match(input_data)
    linearity_coefficient_4 = calibration_variables.linearity_coefficient_4.match(input_data)

    # initiate counter

    counter = PipelineCounter('Linearity', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # check if the zero-read has been identified

        try:

            # get the zero-read and the zero-read error-array from the fits file

            zero_read.from_fits(fits)
            zero_read_error.from_fits(fits)

        except KeyError:

            # locate and rename the zero-read and the zero-read error-array
            # (first non-destructive read - i.e. last in the list of scientific frames)

            fits[sci(fits)[-1]].name = zero_read.keyword
            fits[err(fits)[-1]].name = zero_read_error.keyword

            # get the zero-read and the zero-read error-array from the fits file

            zero_read.from_fits(fits)
            zero_read_error.from_fits(fits)

        # calculate (and save in the data file) the zero-read flux

        zero_read_flux.to_fits(fits, value=zero_read.value - super_zero_read)

        # calculate (and save in the data file) the zero-read flux error-array

        zero_read_flux_error.to_fits(fits, value=np.sqrt(np.abs((ccd_read_noise * ccd_read_noise +
                                                                 ccd_gain * np.abs(zero_read_flux.value))) /
                                                         (ccd_gain * ccd_gain) +
                                                         super_zero_read_error * super_zero_read_error))

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            # just maths, to make the calculation shorter and faster

            c1 = linearity_coefficient_1
            c2 = linearity_coefficient_2
            c3 = linearity_coefficient_3
            c4 = linearity_coefficient_4

            f_r = science
            f_r2 = f_r * f_r

            f_z = zero_read_flux.value
            f_z2 = f_z * f_z
            f_z3 = f_z2 * f_z

            f_rpz = f_r + f_z
            f_rpz2 = f_rpz * f_rpz
            f_rpz3 = f_rpz2 * f_rpz

            f_rmz = f_r * f_z

            f_rzz = f_rpz + f_z

            f_r_err = error
            f_r_err2 = f_r_err * f_r_err

            f_z_err = zero_read_flux_error.value
            f_z_err2 = f_z_err * f_z_err

            # apply the non-linearity correction function

            science = (
                f_rpz * (1.0 + c1 + c2 * f_rpz + c3 * f_rpz2 + c4 * f_rpz3) -
                f_z * (1.0 + c1 + c2 * f_z + c3 * f_z2 + c4 * f_z3)
            )

            # propagate the uncertainties in the error-array

            error = np.sqrt(
                f_r_err2 * ((1.0 + c1 + f_rpz * (2.0 * c2 + f_rpz * (3.0 * c3 + 4.0 * c4 * f_rpz))) ** 2) -
                f_z_err2 * f_r2 * ((2.0 * c2 + 3.0 * c3 * f_rzz + 4.0 * c4 * (f_r2 + 3.0 * f_rmz + 3.0 * f_z2)) ** 2)
            )

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def dark(input_data, splitting=False):
    """
    Dark current correction.

    1) Corrects each NDR for the dark current.



    by subtracting from it the suitable dark current frame
    included in the suitable dark calibration file.
    This frame is selected in order to have the same sub-array, sampling sequence and sample number.

    2) Propagates the uncertainties in the error-array of each NDR.

    The existing error-array is combined with the error-array of the selected dark current (DC_e).

    .. math:: NDR_e \rightarrow \sqrt{NDR_e ^ 2 + DC_e ^ 2}

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for dark current

    """

    # load pipeline and calibration variables to be used

    super_dark = calibration_variables.super_dark.match(input_data)

    # initiate counter

    counter = PipelineCounter('Dark', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            if splitting:
                dark_current = np.array(super_dark[i].data - super_dark[i + 5].data, dtype=float)
                dark_current_error = np.sqrt(super_dark[j].data ** 2 + super_dark[j + 5].data ** 2, dtype=float)
            else:
                dark_current = np.array(super_dark[i].data, dtype=float)
                dark_current_error = np.array(super_dark[j].data, dtype=float)

            # subtract the dark current

            science -= dark_current

            # propagate the uncertainties in the error-array

            error = np.sqrt(error * error + dark_current_error * dark_current_error)

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def gain(input_data):
    """
    Gain variations correction and DN to electrons conversion.

    1) Corrects each NDR for the gain variations among the quarters.

    The WFC3/IR detector consists of four different parts (quadrants) which have different A/D converters. To correct
    for the differences in their gain each NDR is divided by the detector gain pattern (GP), which is included in the
    WFC3 pfl calibration file, and gives the ratio between the gain of each pixel and the average gain of the whole
    detector.

    .. math:: NDR \rightarrow NDR / GP
       
    2) Converts each NDR from DN units to electron units.

    For this correction, each NDR is multiplied by the mean gain of the whole detector (MG), which is included in the
    WFC3 ccd calibration file, and is the average electrons / DN rate of the whole detector.

    .. math:: NDR \rightarrow NDR * MG

    3) Converts the error-array of each NDR from DN units to electron units.

    The same conversion that is applied to each NDR is applied to the each NDR error-array (NDR_e) as well.

    .. math:: NDR_e \rightarrow NDR_e * MG / GP

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for gain variations and converted from DN to electrons

    """

    # load pipeline and calibration variables to be used

    gain_profile = calibration_variables.gain_profile.match(input_data)
    mean_gain = calibration_variables.mean_gain.match(input_data)

    # initiate counter

    counter = PipelineCounter('Gain', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            # correct the NDR for the gain variations

            science /= gain_profile

            # convert the NDR from DN units to electron units

            science *= mean_gain

            # convert the error-array of the NDR from DN units to electron units

            error *= mean_gain / gain_profile

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def sky(input_data, sky_detection_limit=None, splitting=False):
    """
    Sky background correction.

    Corrects each read for the sky background, by subtracting a scaled version of the master sky frame
    included in the suitable sky calibration file.

    The sky calibration file is selected in order to have the same grism (G102 or G141). The scaling factor
    (or sky ratio) is the median of the division between the read and the master sky frame,
    as calculated from a non-illuminated area of the detector.

    The non-illuminated pixels are found based on all the possible differential reads.
    A non-illuminated pixel has to be within mean +/- x * sigma (x = sky_detection_limit)
    in all the possible differential reads, after an eight-pixels moving average soothing
    (this is not affecting the science frame, it is performed only for the detection of the non-illuminated pixels).

    Also propagates the uncertainties in the error arrays.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    sky_detection_limit : float
        detection threshold for the non-illuminated pixels

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for sky background

    """

    # load pipeline and calibration variables to be used

    sky_detection_limit = pipeline_variables.sky_detection_limit.custom(sky_detection_limit)

    sky_background_level = pipeline_variables.sky_background_level.custom()
    sky_frame = pipeline_variables.sky_area.custom()

    master_sky = calibration_variables.master_sky.match(input_data)

    # filter the zero values in the master-sky

    master_sky = np.where(master_sky == 0, 1, master_sky)

    # initiate counter

    counter = PipelineCounter('Sky', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        if not splitting:

            # calculate all the differential combinations of the NDRs (the first read is included as it is)

            differential_science = [fits[sci(fits)[-1]].data]

            for i in range(len(sci(fits))):
                for j in range(i + 1, len(sci(fits))):
                    differential_science.append(np.array(fits[sci(fits)[i]].data - fits[sci(fits)[j]].data))

            # apply an eight-pixel moving average smoothing on the differential combinations

            for i in range(len(differential_science)):
                frame = differential_science[i]
                differential_science[i] = \
                    np.mean([frame,
                             np.roll(frame, 1, 0), np.roll(frame, -1, 0), np.roll(frame, 1, 1), np.roll(frame, -1, 1),
                             np.roll(np.roll(frame, 1, 0), 1, 1), np.roll(np.roll(frame, 1, 0), -1, 1),
                             np.roll(np.roll(frame, -1, 0), 1, 1), np.roll(np.roll(frame, -1, 0), -1, 1)
                             ], 0)

            # detect the sky area in each differential combination

            for i in range(len(differential_science)):

                # calculate the distribution of the flux and fit a gaussian, do not consider pixels at the edges

                frame = differential_science[i][10:-10, 10:-10]
                frame_mean, frame_std = tools.distribution(frame.flatten(), xstep=10.0)[-1][-2:]

                # find the non-illuminated pixels and save them

                differential_science[i] = (np.abs(frame - frame_mean) < sky_detection_limit.value * frame_std)

            # detect the final sky area from all the differential combination

            sky_area = np.where(np.all(np.array(differential_science), 0))
            sky_area = (sky_area[0] + 10, sky_area[1] + 10)

            # save the sky area in the fits file (0 is where the pixels is illuminated only by the sky background)

            sky_map_array = np.ones_like(fits[1].data)
            sky_map_array[sky_area] = 0
            sky_frame.set(sky_map_array)
            sky_frame.to_fits(fits)

            # save the sky detection limit in the fits file

            sky_detection_limit.to_fits(fits)

        else:

            sky_frame.from_fits(fits)
            sky_area = np.where(sky_frame.value == 0)

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            # calculate the sky background level

            sky_background_level.set(np.median((science / master_sky)[sky_area]))

            # subtract the master-sky scaled to the sky background level

            science -= master_sky * sky_background_level.value

            # propagate the uncertainties in the error-array

            # TODO: Consider using the standard deviation of the sky ratio instead of the sky photon noise

            error = np.sqrt(error * error + np.abs(master_sky * sky_background_level.value))

            # save the sky background level in the fits file

            sky_background_level.to_fits(fits, position=i)

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def flat(input_data):
    """
    Flat-field correction.

    1) Corrects each NDR for the flat-field.

    by dividing it by the wavelength-dependent flat-field calculated ... .

    2) Propagates the uncertainties in the error-arrays.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either a single image or a complete data set

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for gain variations

    """

    # load pipeline and calibration variables to be used

    normalised_wavelength_frame = pipeline_variables.normalised_wavelength_frame

    flat_field_coefficient_1 = calibration_variables.flat_field_coefficient_1.match(input_data)
    flat_field_coefficient_2 = calibration_variables.flat_field_coefficient_2.match(input_data)
    flat_field_coefficient_3 = calibration_variables.flat_field_coefficient_3.match(input_data)
    flat_field_coefficient_4 = calibration_variables.flat_field_coefficient_4.match(input_data)

    # initiate counter

    counter = PipelineCounter('Flat', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        # get the wavelength of each pixel from the fits file

        normalised_wavelength_frame.from_fits(fits)

        # calculate the wavelength-dependent flat-field

        flat_field = (flat_field_coefficient_1 +
                      flat_field_coefficient_2 * normalised_wavelength_frame.value +
                      flat_field_coefficient_3 * (normalised_wavelength_frame.value ** 2) +
                      flat_field_coefficient_4 * (normalised_wavelength_frame.value ** 3))

        # exclude outliers

        flat_field = np.where(flat_field < 0.01, 1, flat_field)
        flat_field = np.where(flat_field > 1.1, 1, flat_field)

        # iterate over the NDRs

        for i, j in sci_err(fits):

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            # divide by the flat field

            science /= flat_field

            # propagate the uncertainties in the error-array

            error /= flat_field

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    # return the updated version of the input data

    return input_data


def bpcr(input_data, cr_detection_limit=None, cr_neighbours=None, use_bpcr_fast_mode=None, splitting=False):
    """
    Bad-pixels and cosmic-rays correction.

    Corrects each read for the bad pixels and the cosmic rays.
    The bad pixels are given by the pixels calibration file while the cosmic rays are detected.

    For each pixel, two flags are calculated;
    the median absolute difference from the N (cr_neighbours) horizontally neighbouring pixels (x-flag) and
    the median absolute difference from the N (cr_neighbours) vertically neighbouring pixels (y-flag).
    Given the cosmic rays detection threshold c (cr_detection_limit),
    if a pixel's x-flag is c times larger than the median x-flag in the column (N closest pixels in the column) and
    it's y-flag is c times larger than the median y-flag in the row (N closest pixels in the row),
    it is identified as frame cosmic ray.

    Parameters
    ----------
    input_data : HDUList or DataSet object
        input data, either frame single image or frame complete data set

    cr_detection_limit : float
        detection threshold for cosmic rays

    cr_neighbours : integer
        number of neighbouring pixels from which the x and y flags are calculated

    fast_mode : bool
        whether to avoid correcting all the read or not

    Returns
    -------
    input_data : HDUList or DataSet object
        updated version of the input_data object, corrected for bad pixels and cosmic rays

    """

    # load pipeline and calibration variables to be used

    cr_neighbours = pipeline_variables.cr_neighbours.custom(cr_neighbours)
    cr_detection_limit = pipeline_variables.cr_detection_limit.custom(cr_detection_limit)
    use_bpcr_fast_mode = pipeline_variables.use_bpcr_fast_mode.custom(use_bpcr_fast_mode)

    bpcr_map = pipeline_variables.bpcr_map.custom()

    bad_pixels = calibration_variables.bad_pixels.match(input_data)

    # initiate counter

    counter = PipelineCounter('BPs & CRs', fits_list_size(input_data))

    # iterate over the list of HDUList objects included in the input data

    for fits in fits_list(input_data):

        if not splitting:

            frame = np.array(fits[1].data, dtype=float)

            # calculate the x and y flags from the median absolute difference between each pixel
            # and each N neighbours (N = cr_neighbours)

            x_flag = []
            y_flag = []

            for i in range(1, 1 + int(cr_neighbours.value) / 2):
                x_flag.append(frame - np.roll(frame, i, 1))
                x_flag.append(frame - np.roll(frame, -i, 1))
                y_flag.append(frame - np.roll(frame, i, 0))
                y_flag.append(frame - np.roll(frame, -i, 0))

            x_flag = np.median(x_flag, 0)
            y_flag = np.median(y_flag, 0)

            # calculate the detection limit as me median x_flag per column (from the N closest pixels in the column)
            # and the median y-flag per row (from the N closest pixels in the row)

            x_detection_limit = []
            y_detection_limit = []

            for i in range(1, 1 + int(cr_neighbours.value) / 2):
                x_detection_limit.append(np.roll(x_flag, i, 1))
                x_detection_limit.append(np.roll(x_flag, -i, 1))
                y_detection_limit.append(np.roll(y_flag, i, 0))
                y_detection_limit.append(np.roll(y_flag, -i, 0))

            x_detection_limit = cr_detection_limit.value * np.median(x_detection_limit, 0)
            y_detection_limit = cr_detection_limit.value * np.median(y_detection_limit, 0)

            # characterise as cosmic rays only the pixes which have both their flags above the respective detection limit

            cr_test = np.where((x_flag > x_detection_limit) & (y_flag > y_detection_limit), 1, 0)
            cr_test = np.where(cr_test == 1)

            # local_x_flag = []
            # local_y_flag = []
            #
            # for i in range(1, 1 + int(cr_neighbours.value) / 2):
            #     local_x_flag.append(np.roll(x_flag, i, 0))
            #     local_x_flag.append(np.roll(x_flag, -i, 0))
            #     local_y_flag.append(np.roll(y_flag, i, 1))
            #     local_y_flag.append(np.roll(y_flag, -i, 1))
            #
            # median_x_flag = np.median(local_x_flag, 0)
            # median_y_flag = np.median(local_y_flag, 0)
            #
            # med_x_flag = np.median(np.abs(local_x_flag - median_x_flag), 0)
            # med_y_flag = np.median(np.abs(local_y_flag - median_y_flag), 0)
            #
            # x_detection_limit = cr_detection_limit.value * med_x_flag
            # y_detection_limit = cr_detection_limit.value * med_y_flag
            #
            # # characterise as cosmic rays only the pixes which have both their flags above the respective detection limit
            #
            # cr_test = np.where((np.abs(x_flag - median_x_flag) > x_detection_limit) & (np.abs(y_flag - median_y_flag) > y_detection_limit), 1, 0)
            # cr_test = np.where(cr_test == 1)

            # in the the bad pixels table the bad pixels are characterised by the number -1
            # characterise the cosmic rays by the number 1
            # and save the final bad pixels and cosmic rays map it in the fts file

            bpcr_map_array = np.ones_like(bad_pixels) * bad_pixels
            bpcr_map_array[cr_test] = 1
            bpcr_map.set(bpcr_map_array)
            bpcr_map.to_fits(fits)

            # save process information in the fts file

            cr_neighbours.to_fits(fits)
            cr_detection_limit.to_fits(fits)

        else:

            bpcr_map.from_fits(fits)

        # scan or no scan

        data = np.sum(fits[1].data, 0)
        model = tools.box(np.arange(len(data)), len(data) / 2, np.max(data), 45., 40.)
        dx = np.argmax(np.convolve(data, model)) - np.argmax(np.convolve(model, model))
        x_lim1, x_lim2 = len(data) / 2 + dx - 55, len(data) / 2 + dx + 55

        # detect the approximate vertical position of the final spectrum

        data = fits[1].data
        y1 = 5 + np.median(np.argmax(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
        y2 = np.median(np.argmin(data[5:, x_lim1:x_lim2] - data[:-5, x_lim1:x_lim2], 0))
        if abs(y2 - y1) <= 1:
            y1 = 2 + np.median(np.argmax(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
            y2 = np.median(np.argmin(data[2:, x_lim1:x_lim2] - data[:-2, x_lim1:x_lim2], 0))
        final_y_lim1, final_y_lim2 = np.sort([int(round(y1)), int(round(y2))])

        if abs(final_y_lim1 - final_y_lim2) > 1:
            final_scan = True
        else:
            final_scan = False

        # correct each read for the bad pixels and the cosmic rays

        if len(sci_err(fits)) == 1:
            corr = [sci_err(fits)[0]]
        else:
            if use_bpcr_fast_mode.value:
                corr = [sci_err(fits)[0], sci_err(fits)[-1]]
            else:
                corr = sci_err(fits)

        for i, j in corr:

            # get the NDR and its error-array from the data file

            science = np.array(fits[i].data, dtype=float)
            science_old = np.array(fits[i].data, dtype=float)
            error = np.array(fits[j].data, dtype=float)

            if final_scan:

                # find the clean pixels location and values

                grid_x, grid_y = np.meshgrid(np.arange(len(science)), np.arange(len(science[0])))
                clean_pixs = np.where(bpcr_map.value == 0)
                points = np.swapaxes(np.roll(clean_pixs, 1, 0), 0, 1)
                values = fits[i].data[clean_pixs].flatten()

                # replace the frame with the interpolated values
                # after this process the clean pixels remain the same
                # and the bad pixels and cosmic rays are replaced

                # TODO: this function is slow, try and write a faster process using np.arrays

                science = griddata(points, values, (grid_x, grid_y), method='cubic')

            else:
                bpcr_map.value[:, 0] = 0
                bpcr_map.value[:, -1] = 0
                for fits_line in range(len(fits[i].data)):
                    clean_pixs = np.where(bpcr_map.value[fits_line] == 0)
                    points = clean_pixs[0]
                    values = fits[i].data[fits_line][clean_pixs]
                    line_model = interp1d(points, values)
                    science[fits_line] = line_model(np.arange(len(science[fits_line])))

            # propagate the uncertainties in the error array assuming that the difference between
            # the original the final values is and additive correction

            # TODO: This calculation is not the most optimal. A better approach would be to assume that the
            # TODO: interpolated value is a new measurement with an uncertainty originated from the interpolation.
            # TODO: The photon noise (square root of the interpolated value) should be added to that uncertainty to
            # TODO: calculate the final uncertainty.

            error = np.sqrt(error ** 2 + np.abs(science - science_old))

            # update the NDR and its error-array in the data file

            fits[i].data = science
            fits[j].data = error

        # update counter

        counter.update()

    return input_data
