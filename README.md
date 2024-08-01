# RGB2SCA
These MATLAB scripts allow for snow cover maps to be produced interactively from optical (RGB) imagery. It was designed to work with RGB orthomosaics produced using Uncrewed Aerial Systems (UAS or drones).

Snow has high reflectance across the full visible portion of the electromagnetic spectrum, but it is highest in blue wavelengths (450 – 490 nm) and comparatively lowest in red wavelengths (590 – 650 nm, Dozier 1989). Because of this, blue band thresholding has been effective for delineating SCA (Thaler et al., 2023) as well as band differences between color channels. Testing of various band ratios, band differences, and single-band thresholds at multiple New Hampshire, U.S.A study sites concluded that the normalized difference between the blue and red bands (blue – red, BRd) was shown to most effectively delineate both shaded and unshaded snow. The magnitude of pixel intensity was also considered to ensure that bright white regions, like snow, were also captured.

##Binary SCA maps are produced by optimizing thresholds of combined pixel intensity and the difference between the blue and red color channels. The processing workflow to achieve this is as follows:

(1) User inputs path to an image (e.g., RGB orthomosaic)

(2) The image is loaded, and the user selects a small reference area and manually creates a snow vs. no snow reference

(3) Step 2 is repeated until the desired number of reference areas are produced

(4) The mean F1-score across a randomly selected subset (~75%) of these reference areas is then maximized by tuning a weighting metric (0 to 1) and snow cover metric (SCM) threshold (-0.1 to 1) and comparing the reference map to the produced map

(5) Performance statistics and rasters are automatically produced and stored to the same path as the input image

Further documentation and information on these scripts is provided as code comments.
