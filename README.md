![PyPI Downloads](https://static.pepy.tech/badge/iraclis)

# Iraclis

<img src="https://github.com/ucl-exoplanets/Iraclis/blob/master/logo.jpg" width="25%">

Analysis pipeline for **HST/WFC3** spectroscopic observations of exoplanet transits and eclipses

A complete package to analyse single-object spatially scanned spectroscopic observations of extrasolar planets 
obtained with the near-infrared grisms (G102, G141) of the **Wide Field Camera 3** on-board the 
**Hubble Space Telescope**. 

Includes:
* Reduction of the raw frames.
* Calibration of the position of the total spectrum and the different spectral elements.
* Extraction of the total flux and the flux per spectral element.
* Fitting of the white and the spectral light-curves.

Currently, fitting can be applied only on single-visit light-curves but in the next version it will
be updated to fit also multiple visits of the same target simultaneously.


## References

* Tsiaras et al. 2016a, [A New Approach to Analyzing HST Spatial Scans: The Transmission Spectrum of HD 209458 b](http://iopscience.iop.org/article/10.3847/0004-637X/832/2/202), ApJ, 832, 202. 
* Tsiaras et al. 2016b, [Detection of an Atmosphere Around the Super-Earth 55 Cancri e](http://iopscience.iop.org/article/10.3847/0004-637X/820/2/99), ApJ, 820, 99.
* Tsiaras et al. 2018, [A Population Study of Gaseous Exoplanets](http://iopscience.iop.org/article/10.3847/1538-3881/aaaf75), AJ, 155, 156.


## License

This work is licensed under the Creative Commons Attribution 4.0 International License.

Copyright (c) 2018-2021 Angelos Tsiaras

Please pay attention to Section 3 of the license and do:
- retain identification of the creators by including the above listed references in future work and publications,
- indicate if You modified the Licensed Material and retain an indication of any previous modifications.

To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ 
or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.


## Installation

Open a terminal and type: `pip install iraclis`  


## Updates

#### v.1.5.4 - Improvements and bugs fixes
- Replace np.product with np.prod 
- Installation issues fixed
- Improved detection of bright and faint stars in the direct image
- Improved eclipse mid-time prediction
- Saving figures to PDF instead of EPS
- Bug fixed: version 1.5 introduced a bug that was muting wavelength-dependent LDCs

#### v1.5.0 - Updates related to the new PyLightcurve v4
- Automated detection of planet parameters
- Conversion of time of BJD_TDB,instead of HJD_UTC


## Usage

#### Getting HST/WFC3 data

Each transit or eclipse dataset, obtained with WFC3, includes many spectroscopic images (using one of the two 
spectroscopic grisms: G102 or G141) and also, at least, one direct (undispersed) image (using one of the 15 WFC3 imaging 
filters: F098W, F132N, F140W, F126N, F153M, F167N, F139M, F164N, F127M, F160W, F128N, F125W, F130N, F110W, F105W).

To use Iraclis, download all the spectroscopic data in the RAW format and the imaging data in the FLT format
from the [MAST archive](https://archive.stsci.edu/hst/search.php). Keep all your files together in one directory.

#### Analysing HST/WFC3 data

Within a bash terminal:
- you can process a dataset using a parameters file by typing: `iraclis -p path_to_parameters_fle`
- or using default values for the parameters by typing: `iraclis -d path_to_data_directory`

Within a python terminal, first import the package: `import iraclis`
- you can process a dataset using a parameters file by typing: `iraclis.process_visit(parameters_file='path_to_parameters_fle')`
- or using default values for the parameters by typing: `iraclis.process_visit(data_directory='path_to_data_directory')`

#### Setting up a parameters file

A parameters file is a simple txt file that contains all the important paramets controlling the data anlaysis process 
of a WFC3 dataset. Below, you can find the description of such a file (this file is included in the examples).

`data_directory` **`iraclis_test_dataset_hatp26b_raw_data`**  
(diectory path) The path of the directory where the dataset is stored.

`output_directory_copy` **`False`**  
(diectory name/False) If you wish to copy your results in a new directory, give here its name. The default results are 
stored in the “results” directory. This is useful if you wish to analyse the same dataset multiple times with different 
fitting parameters.

`reduction` **`True`**  
(True/False) Whether to reduce the data or not. You can set it to False if you wish to re-run the analysis starting 
from a later stage.

`splitting` **`True`**  
(True/False) Whether to split the data or not. In this case the reduced data will be splitted into the differential 
frames N - N-1, rather than final-initial.

`extraction` **`False`**  
(True/False) Whether to extract the light-curves or not. You can set it to False if you wish to re-run the analysis 
starting from a later stage.

`splitting_extraction` **`True`**  
(True/False) Whether to extract the light-curves from the splitted data or not. You can set it to False if you wish to 
re-run the analysis starting from a later stage. Note: set only one of the extraction or splitting_extraction to true. 
If you wish to have both, re-run the analysis with different values for these parameters.

`fitting_white` **`True`**  
(True/False) Whether to fit the white light-curve or not. You can set it to False if you wish to re-run the analysis 
starting from a later stage.

`fitting_spectrum` **`True`**  
(True/False) Whether to fit the spectral light-curves or not. 

`target_x_offset` **`0`**  
`target_y_offset` **`0`**  
(number) In the rare case that in the field of view there is a star brighter that your target, give here the difference 
in their positions. This is important because the code automatically identifies the brightest star as the target 
(in 99% of the observations, this is the case).

`aperture_lower_extend` **`-20`**  
(number) Vertical starting position of the extraction box. Use negative value. -20 means that the extraction box will 
start 20 pixels below the spectrum.

`aperture_upper_extend` **`20`**  
(number) Vertical final position of the extraction box. Use positive value. 20 means that the extraction box will stop 
20 pixels above the spectrum.

`extraction_method` **`gauss`**  
(gauss/integral) There are two available extraction methods: gauss (where the flux extracted is the convolution of the 
spectrum with a gaussian) or integral (where the flux extracted is calculated as the integral of the spectrum inside the 
extraction aperture)

`extraction_gauss_sigma` **`47`**  
(number) Useful only in the case of gauss extraction. This is the sigma of the gaussian used, in Angstrom. The default number of 
47 is approximately equal to one pixel.

`method` **`claret`**  
(claret/linear/quad/sqrt) Limb-darkening method to be used. Choose between: claret, linear, quad, sqrt.

`white_lower_wavelength` **`default`**  
`white_upper_wavelength` **`default`**  
(number/default) Right and left edges of the extracted white light curve, in Angstroms. 

`white_ldc1` **`default`**  
(number/default) First limb-darkening coefficients for the white light-curve.

`white_ldc2` **`default`**  
(number/default) Second limb-darkening coefficients for the white light-curve. Will not be used if the linear method is 
chosen.

`white_ldc3` **`default`**  
(number/default) Third limb-darkening coefficients for the white light-curve. Will not be used if the linear, the quad 
or the sqrt method is chosen.

`white_ldc4` **`default`**  
(number/default) Forth limb-darkening coefficients for the white light-curve. Will not be used if the linear, the quad 
or the sqrt method is chosen.

**Comment**: When the above parameters are set to default, the white light-curve limits will be 10880.0 - 16800.0 
Angstroms for G141 and 8000 - 11500 Angstroms for G102.

`bins_file` **`default_high`**  
(file path/default_high/default_low/default_vlow) Path to the bins file.

**Comment**: A bins file should contain 4 columns: 1 - bins lower edge in Angstrom, 2 - bins upper edge in Angstrom, 
3 - first limb darkening coefficient, 4 - second limb darkening coefficient, 5 - third limb darkening coefficient, 
6 - forth limb darkening coefficient. An example of a bins file can be found in iraclis_test_dataset_hatp26b_bins.txt

`planet` **`auto`**  
(name - no spaces/auto) Planet name.

`star_teff` **`auto`**  
(number/auto) Stellar temperature, used if the limb-darkening coefficients are set to auto.

`star_logg` **`auto`**  
(number/auto) Stellar log(g), used if the limb-darkening coefficients are set to auto.

`star_meta` **`auto`**  
(number/auto) Stellar metallicity, used if the limb-darkening coefficients are set to auto.

`rp_over_rs` **`auto`**  
(number/auto) Planet-to-star radius ratio. This parameters is always fitted both for the white 
and the spectral light-curves in cases of transits.

`fp_over_fs` **`auto`**  
(number/auto) Planet-to-star flux ratio. This parameters is always fitted both for the white and 
the spectral light-curves in cases of eclipses.

`period` **`auto`**  
(number/auto) Period of the planetary orbit in days. Always fixed.

`sma_over_rs` **`auto`**  
(number/auto) Semi-major axis of the planetary orbit.

`eccentricity` **`auto`**  
(number/auto) Eccentricity of the planetary orbit. Always fixed.

`inclination` **`auto`**  
(number/auto) Inclination of the planetary orbit, in degrees.

`periastron` **`auto`**  
(number/auto) Periastron of the planetary orbit in degrees. Always fixed.

`mid_time` **`auto`**  
(number/auto) Mid-transit time of the planetary orbit, in BJD_TDB.
This should be the mid-transit time, even if the observation is an eclipse one.

**Comment**: You can set any of the above 12 parameters to auto, to use the data from a built-in catalogue.

`apply_up_down_stream_correction` **`False`**  
(True/False) Whether to correct for the up-stream/down-stream effect or not. Useful only in cases of fast scans that 
cross the line between the upper two and lower two quadrants of the detector.

`exclude_initial_orbits` **`1`**  
(number) Number of HST orbits to be removed from the begging of the visit. Usually set to 1.

`exclude_final_orbits` **`0`**  
(number) Number of HST orbits to be removed from the end of the visit. Usually set to 0.

`exclude_initial_orbit_points` **`0`**  
(number) Number of HST exposures to be removed from the begging of each HST-orbit.

`mcmc_iterations` **`300000`**  
(number) Number of emcee iterations for the white light-curve fitting

`mcmc_walkers` **`50`**  
(number) Number of emcee wakers for the white light-curve fitting

`mcmc_burned_iterations` **`200000`**  
(number) Number of emcee burned iterations for the white light-curve fitting

`spectral_mcmc_iterations` **`60000`**  
(number) Number of emcee iterations for the spectral light-curve fitting

`spectral_mcmc_walkers` **`50`**  
(number) Number of emcee walkers for the spectral light-curve fitting

`spectral_mcmc_burned_iterations` **`10000`**  
(number) Number of emcee burned iterations for the spectral light-curve fitting

`first_orbit_ramp` **`True`**  
(True/False) Whether to fit for different ramp coefficients for the first HST orbit in the analysis (after excluding 
exclude_initial_orbits orbits) or not.

`second_order_ramp` **`False`**  
(True/False) Whether to fit for a quadratic visit-long ramp or not.

`mid_orbit_ramps` **`True`**  
(True/False) Whether to fit for mid-orbit ramps caused by buffer-dumps or not.

`fit_ldc1` **`False`**  
(True/False) Whether to fit for the first limb-darkening coefficient or not. The same will be applied both for the 
white and the spectral light-curves.

`fit_ldc2` **`False`**  
(True/False) Whether to fit for the second limb-darkening coefficient or not. The same will be applied both for 
the white and the spectral light-curves. Will not be used if the linear method is chosen.

`fit_ldc3` **`False`**  
(True/False) Whether to fit for the third limb-darkening coefficient or not. The same will be applied both for 
the white and the spectral light-curves. Will not be used if the linear, the quad or the sqrt method is chosen.

`fit_ldc4` **`False`**  
(True/False) Whether to fit for the forth limb-darkening coefficient or not. The same will be applied both for 
the white and the spectral light-curves. Will not be used if the linear, the quad or the sqrt method is chosen.

`fit_sma_over_rs` **`False`**  
(True/False) Whether to fit for the semi-major axis of the planetary orbit or not. This is fitted only on the white 
light-curve.

`fit_inclination` **`False`**  
(True/False) Whether to fit for the inclination of the planetary orbit or not. This is fitted only on the white 
light-curve.

`fit_mid_time` **`True`**  
(True/False) Whether to fit for the mid-transit or mid-eclipse time of the planetary orbit or not. This is fitted only on the white 
light-curve.

#### Setting up a bins file

A bins file is a simple txt file that contains three to six columns, indicating for each spectral light curve: 
a-b. Right and left edges of the extracted spectral light curve, in Angstroms, c-f. limb-darkening coefficients for 
the spectral light-curve. Below, you can find the description of such a file (this file is included in the examples).

`11108 11416 0.985047 -1.385670 1.781030 -0.7267230`  
`11416 11709 0.949097 -1.266470 1.640630 -0.6754740`  
`11709 11988 0.928715 -1.195690 1.553730 -0.6452910`  
`11988 12257 0.903069 -1.109910 1.456180 -0.6107730`  
`12257 12522 0.878225 -1.023230 1.361070 -0.5780620`  
`12522 12791 0.859841 -0.950740 1.274760 -0.5481460`  
`12791 13058 0.849884 -0.896126 1.203900 -0.5267150`  
`13058 13321 0.832077 -0.833290 1.125660 -0.4941230`  
`13321 13586 0.809188 -0.726211 0.991314 -0.4438480`  
`13586 13860 0.795028 -0.641081 0.872551 -0.3971340`  
`13860 14140 0.788556 -0.586294 0.802106 -0.3739860`  
`14140 14425 0.784454 -0.561833 0.775730 -0.3685690`  
`14425 14719 0.772069 -0.460859 0.627183 -0.3091360`  
`14719 15027 0.772703 -0.404730 0.517165 -0.2597780`  
`15027 15345 0.772846 -0.296638 0.327104 -0.1754390`  
`15345 15682 0.798113 -0.256525 0.198611 -0.1108030`  
`15682 16042 0.848830 -0.376905 0.274511 -0.1251750`  
`16042 16432 0.894871 -0.410984 0.233093 -0.0939863`

#### Testing the code

In the examples you will find a short python script that downloads a test dataset from the MAST 
archive. There, you can also find the parameters file and the bins file described above. The test dataset is a transit 
of HATP-26 b from the observing proposal 14260 (PI: Deming Drake). 


## Products

The final product is a pickle file named "fitting_results.pickle" and can be found in the "results" directory 
(or the output_directory_copy that you have set in the parameters file. To open the pickle file you will need python 
installed in your computer and the packages pickle and numpy.

The command to load a pickle file is: `database = pickle.load(open(‘database_file.pickle’))`

This dictionary contains the three main dictionaries that include all the results for a particular dataset. These are:

`database['lightcurves']` extracted light-curves and transit fitting data  
`database['spectrum']` planetary spectrum  

The `database['lightcurves']` dictionary contains all the information on the extracted white and spectral 
light-curves and their fitting, and has the following content:

`database['lightcurves']['white']` white light curve and fitting  
`database['lightcurves']['bin_01']` first wavelength channel  
`database['lightcurves']['bin_02']` second wavelength channel  
etc... 

The `database['spectra']` dictionary contains the final extracted planetary spectrum, and has the following structure:

`database['spectrum']['wavelength']` mean wavelength of the different channels (μm)  
`database['spectrum']['width']` wavelength width of the different channels (μm)  
`database['spectrum']['depth']` transit depth in each wavelength channel  
`database['spectrum']['error']` uncertainty in the transit depth in each wavelength channel  

#### White light curve and fitting

The `database['lightcurves']['white']` dictionary has the following structure:

`database['lightcurves']['white']['limb_darkening']['method']` limb darkening method used  
`database['lightcurves']['white']['wavelength']['lambda1']` lower wavelength limit of the band in Angstroms  
`database['lightcurves']['white']['wavelength']['lambda2']` upper wavelength limit of the band in Angstroms  
`database['lightcurves']['white']['wavelength']['lambda_mean']` mean wavelength of the band in Angstroms  
`database['lightcurves']['white']['wavelength']['lambda_width']` width of the band in Angstroms  
`database['lightcurves']['white']['exposure']['exp_time']` exposure time for each data point in seconds  
`database['lightcurves']['white']['exposure']['model_resolution']` sub-exposure time used to model each data point in seconds  
`database['lightcurves']['white']['input_time_series']['x_shift']` horizontal shifts during the observation  
`database['lightcurves']['white']['input_time_series']['x_shift_error']` uncertainty in the horizontal shifts during the observation  
`database['lightcurves']['white']['input_time_series']['y_shift']` vertical shifts during the observation  
`database['lightcurves']['white']['input_time_series']['y_shift_error']` uncertainty in the horizontal shifts during the observation  
`database['lightcurves']['white']['input_time_series']['sky']` sky background level during the observation  
`database['lightcurves']['white']['input_time_series']['scan']` scan direction (1 for forward scans, -1 for reverse scans)  
`database['lightcurves']['white']['input_time_series']['hjd']` heliocentric Julian date during the observation  
`database['lightcurves']['white']['input_time_series']['raw_lc']` raw white light-curve  
`database['lightcurves']['white']['input_time_series']['raw_lc_error']` uncertainty in the raw white light-curve  
`database['lightcurves']['white']['output_time_series']['phase']` orbital phase  
`database['lightcurves']['white']['output_time_series']['systematics']` best-fit model for the systematics  
`database['lightcurves']['white']['output_time_series']['detrended_lc']` de-trended white light-curve  
`database['lightcurves']['white']['output_time_series']['transit']` best-fit model for the transit  
`database['lightcurves']['white']['output_time_series']['residuals']` fitting residuals  
`database['lightcurves']['white']['statistics']['res_std']` standard deviation of the residuals  
`database['lightcurves']['white']['statistics']['res_autocorr']` autocorrelation function of the residuals  
`database['lightcurves']['white']['statistics']['corr_variables']` fitted variables  
`database['lightcurves']['white']['statistics']['corr_matrix']` correlation matric of the fitted variables  
`database['lightcurves']['white']['parameters']['ldc_1']` first limb-darkening coefficient  
`database['lightcurves']['white']['parameters']['ldc_2']` second limb-darkening coefficient  
`database['lightcurves']['white']['parameters']['ldc_3']` third limb-darkening coefficient  
`database['lightcurves']['white']['parameters']['ldc_4']` forth limb-darkening coefficient  
`database['lightcurves']['white']['parameters']['rp']` planetary radius relative to the stellar radius  
`database['lightcurves']['white']['parameters']['fp']` planetary flux relative to the stellar flux (useful only for eclipses)  
`database['lightcurves']['white']['parameters']['P']` orbital period in days  
`database['lightcurves']['white']['parameters']['a']` orbital semi-major axis relative to the stellar radius  
`database['lightcurves']['white']['parameters']['e']` orbital eccentricity  
`database['lightcurves']['white']['parameters']['i']` orbital inclination in degrees  
`database['lightcurves']['white']['parameters']['omega']` orbital argument of periastron in degrees  
`database['lightcurves']['white']['parameters']['t_0']` mit-transit time in HJD  
`database['lightcurves']['white']['parameters']['n_w_for']` normalization factor for the forward scans  
`database['lightcurves']['white']['parameters']['n_w_rev']` normalization factor for the reverse scans  
`database['lightcurves']['white']['parameters']['r_a1']` linear term of the long-term ramp  
`database['lightcurves']['white']['parameters']['r_a2']` quadratic term of the long-term ramp  
`database['lightcurves']['white']['parameters']['r_b1']` amplitude of the exponential short-term ramp  
`database['lightcurves']['white']['parameters']['r_b2']` decay of the exponential short-term ramp  
`database['lightcurves']['white']['parameters']['for_b1']` amplitude of the exponential short-term ramp for the first orbit  
`database['lightcurves']['white']['parameters']['for_b2']` amplitude of the exponential short-term ramp for the first orbit  
`database['lightcurves']['white']['parameters']['mor_b1']` amplitude of the exponential short-term ramp after a buffer dump  
`database['lightcurves']['white']['parameters']['mor_b2']` amplitude of the exponential short-term ramp after a buffer dump  

Each `database['lightcurves']['white']['parameters']['xx']` element includes also the following keys:

`database['lightcurves']['white']['parameters']['xx']['name']` name  
`database['lightcurves']['white']['parameters']['xx']['initial']` initial value (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['min_allowed']` minimum of the prior (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['max_allowed']` maximum of the prior (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['trace']` mcmc trace (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['trace_bins']` mcmc trace bins (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['trace_counts']` mcmc trace distribution (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['value']` final value  
`database['lightcurves']['white']['parameters']['xx']['m_error']` final -error (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['p_error']` final +error (None if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['print_name']` name shown in plots  
`database['lightcurves']['white']['parameters']['xx']['print_value']` value shown in plots  
`database['lightcurves']['white']['parameters']['xx']['print_m_error']` -error shown in plots (- if not fitted)  
`database['lightcurves']['white']['parameters']['xx']['print_p_error']` +error shown in plots (- if not fitted)  

#### Spectral light curves and fitting

Each `database['light_curves']['bin_yy']` dictionary has the following structure:

`database['lightcurves']['bin_yy']['limb_darkening']['method']` limb darkening method used  
`database['lightcurves']['bin_yy']['wavelength']['lambda1']` lower wavelength limit of the band in Angstroms  
`database['lightcurves']['bin_yy']['wavelength']['lambda2']` upper wavelength limit of the band in Angstroms  
`database['lightcurves']['bin_yy']['wavelength']['lambda_mean']` mean wavelength of the band in Angstroms  
`database['lightcurves']['bin_yy']['wavelength']['lambda_width']` width of the band in Angstroms  
`database['lightcurves']['bin_yy']['exposure']['exp_time']` exposure time for each data point in seconds  
`database['lightcurves']['bin_yy']['exposure']['model_resolution']` sub-exposure time used to model each data point in seconds  
`database['lightcurves']['bin_yy']['input_time_series']['x_shift']` horizontal shifts during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['x_shift_error']` uncertainty in the horizontal shifts during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['y_shift']` vertical shifts during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['y_shift_error']` uncertainty in the horizontal shifts during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['sky']` sky background level during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['scan']` scan direction (1 for forward scans, -1 for reverse scans)  
`database['lightcurves']['bin_yy']['input_time_series']['hjd']` heliocentric Julian date during the observation  
`database['lightcurves']['bin_yy']['input_time_series']['raw_lc']` raw spectral light-curve  
`database['lightcurves']['bin_yy']['input_time_series']['white_raw_lc']` raw white light-curve  
`database['lightcurves']['bin_yy']['input_time_series']['relative_lc']` relative spectral light-curve (raw relative spectral light-curve divided by the raw white light-curve)
`database['lightcurves']['bin_yy']['input_time_series']['relative_lc_error']` uncertainty in the relative spectral light-curve
`database['lightcurves']['bin_yy']['input_time_series']['white_model']` transit model for the white light-curve  
`database['lightcurves']['bin_yy']['output_time_series']['phase']` orbital phase  
`database['lightcurves']['bin_yy']['output_time_series']['systematics']` best-fit model for the systematics  
`database['lightcurves']['bin_yy']['output_time_series']['detrended_lc']` de-trended white light-curve  
`database['lightcurves']['bin_yy']['output_time_series']['transit']` best-fit model for the transit  
`database['lightcurves']['bin_yy']['output_time_series']['residuals']` fitting residuals  
`database['lightcurves']['bin_yy']['statistics']['res_std']` standard deviation of the residuals  
`database['lightcurves']['bin_yy']['statistics']['res_autocorr']` autocorrelation function of the residuals  
`database['lightcurves']['bin_yy']['statistics']['corr_variables']` fitted variables  
`database['lightcurves']['bin_yy']['statistics']['corr_matrix']` correlation matric of the fitted variables  
`database['lightcurves']['bin_yy']['parameters']['ldc_1']` first limb-darkening coefficient  
`database['lightcurves']['bin_yy']['parameters']['ldc_2']` second limb-darkening coefficient  
`database['lightcurves']['bin_yy']['parameters']['ldc_3']` third limb-darkening coefficient  
`database['lightcurves']['bin_yy']['parameters']['ldc_4']` forth limb-darkening coefficient  
`database['lightcurves']['bin_yy']['parameters']['rp']` planetary radius relative to the stellar radius  
`database['lightcurves']['bin_yy']['parameters']['fp']` planetary flux relative to the stellar flux (useful only for eclipses)  
`database['lightcurves']['bin_yy']['parameters']['P']` orbital period in days  
`database['lightcurves']['bin_yy']['parameters']['a']` orbital semi-major axis relative to the stellar radius  
`database['lightcurves']['bin_yy']['parameters']['e']` orbital eccentricity  
`database['lightcurves']['bin_yy']['parameters']['i']` orbital inclination in degrees  
`database['lightcurves']['bin_yy']['parameters']['omega']` orbital argument of periastron in degrees  
`database['lightcurves']['bin_yy']['parameters']['t_0']` mit-transit time in HJD  
`database['lightcurves']['bin_yy']['parameters']['n_l_for']` normalization factor for the forward scans  
`database['lightcurves']['bin_yy']['parameters']['n_l_rev']` normalization factor for the reverse scans  
`database['lightcurves']['bin_yy']['parameters']['r_a1']` linear term of the long-term ramp   

Each `database['lightcurves']['bin_yy']['parameters']['zz']` element includes also the following keys:

`database['lightcurves']['bin_yy']['parameters']['zz']['name']` name  
`database['lightcurves']['bin_yy']['parameters']['zz']['initial']` initial value (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['min_allowed']` minimum of the prior (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['max_allowed']` maximum of the prior (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['trace']` mcmc trace (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['trace_bins']` mcmc trace bins (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['trace_counts']` mcmc trace distribution (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['value']` final value  
`database['lightcurves']['bin_yy']['parameters']['zz']['m_error']` final -error (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['p_error']` final +error (None if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['print_name']` name shown in plots  
`database['lightcurves']['bin_yy']['parameters']['zz']['print_value']` value shown in plots  
`database['lightcurves']['bin_yy']['parameters']['zz']['print_m_error']` -error shown in plots (- if not fitted)  
`database['lightcurves']['bin_yy']['parameters']['zz']['print_p_error']` +error shown in plots (- if not fitted)  


## BUGS!!!

For any issues and bugs please send an E-mail at [atsiaras@star.ucl.ac.uk](atsiaras@star.ucl.ac.uk).
