# Iraclis

Analysis pipeline for **HST/WFC3** spectroscopic observations of exoplanet transits and eclipses

A complete package to analyse single-object spatially scanned spectroscopic observations of extrasolar planets 
obtained with the near-infrared grisms (G102, G141) of the **Wide Field Camera 3** on-board the **Hubble Space Telescope**. 

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

Copyright (c) 2018 Angelos Tsiaras

Please pay attention to Section 3 of the license and do:
- retain identification of the creators by including the above listed references in future work and publications,
- indicate if You modified the Licensed Material and retain an indication of any previous modifications.

To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ 
or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.


## Installation

Open a terminal and type: `pip install iraclis`

Or download, unzip, `cd` in the unzipped directory, and type: `python setup.py install`

You will find a test dataset, including instructions on how to use the pipeline here (62.3 MB): 
https://www.dropbox.com/s/l3jcz8dzmw5mwg6/iraclis_test_dataset_hatp26b.zip?dl=0

For any issues please send an E-mail at: atsiaras@star.ucl.ac.uk