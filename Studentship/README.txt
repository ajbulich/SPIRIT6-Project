10/2/23 - Anthony Bulich

Final iteration of code. This README will outline which .py files do what and how to use them. 
Docstrings have also been written for the functions, to describe the parameters they require.

The file structure is that all the data is stored inside SPIRIT6-Observations/GoodDarks/MultipleExposures within this folder, every night of observations has its own folder, yyyymmdd, which contains a separate folder for each object imaged that night

PhotutitilsTutorials folder contains all the code. The folder structure has been left in-tact, with 2 example folders that still have data inside, to help understand the overall structure and file formats.

One thing to note is that all the object names for any inputs into any function are not exactly 
as written in their proper scientific name. If there are any spaces, they are removed, and if there
are any hyphens, only the first chunk before the hyphen is used.
e.g. MCT 0550-4911 becomes MCT0550, WD 0830-535 becomes WD0830. This is only when used in functions
to make everything easier.


The file format for 90secExposures{objectName}_{FilterType}_{date}.txt is just the path to all the 
images, each path on a new line

The file format for {objectName}Stars.txt is - 
RA Dec Magnitude(B) Magnitude(V)

for every star in the field

These two file types must be manually written.


FinalAperturePhotometry.py - this files contains the function annulusPhotometry, which performs 
                             the aperture photometry used in the pipeline.
FluxAccuracy.py - creates a plot for flux accuracy for those objects on the given nights. most
                  importantly however, it writes the files (Mean/Median)CollectiveMagnitudeResults(B/V).txt
                  which are used to create flux accuracy histograms for all dates combined
FluxCalibration1.py - Contains the function FluxCalibration which performs the flux calibration as the 
                      last step of the overall pipeline. Also contains RepeatedFluxCalibration
FluxHistogram.py - contains the function createFluxHistogram, which produces the cumulative flux accuracy
                   histogram for either mean or median stacked data. FluxAccuracy.py must be run first
                   for every date/object combination.
FluxPipeline.py - this file contains the main pipeline that reprojects, stacks, performs photometry and 
                  then flux calibrates the stars in an image. It contains FluxCalibrationAnalysis as well
                  as automatedAnalysis. The latter is most used as it will create the incremented stacked images
                  as well as do both median and mean. 
HistogramComparison.py - Creates histograms of incrementally stacked images for the purpose of determined noise profile
NoiseTimeUpdated.py - creates the noise vs time plots for a given object and set of dates. The main pipeline must
                      be run first however, to create the incrementally stacked images so that the noise vs time
                      plot can be created. This is the job of automatedAnalysis.
projectMultiple.py - contains the code for image reprojection, first step of pipeline.
psfForStackedImages.py - contains SingleImagePSF, which outputs the sigma values for the stars in a single image
                         also contains PSFStackedPhotometry, which performs the PSF photometry on the stacked images
                         as part of the main pipeline.
SeeingHistogram.py - contains addData, which adds the seeing data for specific dates to a file, in mode 1.
                     then mode 2 access the file for every date and appends it to one big file.
                     also contains createHistogram which can create a histogram for a specific date, or for
                     all dates cumulatively.
SeeingPlot.py - creates the seeing over a single night plots. e.g plots FWHM from 9pm to 3am. Function is plotSeeing.
stackMultiple.py - contains the functions that perform the image stacking for the main pipeline, as well as a modified 
                   function that performs incremental stacking.
StackWeight.py - contains StackWeighting, which determines the weighting array for a group of images.


