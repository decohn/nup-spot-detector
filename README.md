# nup-spot-detector

Code to compute the integrated intensity of isolated Nup foci, designed for the Reyes-Lamothe Lab.

## Table of Contents

I. **Project Goals**

II. **Outline of Algorithm**

III. **Using the Program**

	a. Setting up the Configuration File 
	b. Optimizing Parameters
	c. Understanding the Output
	d. Applications of the Output
	
IV. **Other References**

	a. Project Proposal (October 15, 2018)
	b. Lab Meeting Slides (November 15, 2018)
	c. Interim Project Report (December 4, 2018)
	d. Lab Meeting Slides (February 27, 2019)
	e. Final Lab Meeting Slides (April 18, 2019)
	
V. **Authorship and Acknowledgements**

	a. Authorship
	b. Acknowledgements
	
## I. Project Goals

This script was desgined to accomplish the goals of part one of the replisome protein copy number estimation project. Specifically, it will identify the average fluorescence intensity per molecule of mNeonGreen (or any other fluorescent tag) by leveraging the known stoichiometry of the nuclear pore complex. Given a group of image stacks for a particular yeast strain in which a single NPC protein is fluorescently tagged, as well as the manually identified locations of isolated NPCs in those stacks, the script will report the average fluorescence intensity above background of an NPC. In principle, this will also account for autofluorescence. The results from multiple runs of this script on images from different yeast strains can be used to construct a standard curve of NPC fluorescence intensity vs. tagged protein stoichiometry within the NPC, which can be used to estimate the fluorescence intensity of a single fluorescent tag molecule.

## II. Outline of Algorithm

For a visual representation of the below outline, please see the flowchart included in this repository.

1. Read input from the SNSDconfig.txt file, namely the names of the image files being used and the locations of the isolated Nup spots within each image stack.

2. Read the first image stack.

FOR EACH isolated Nup spot:

3a. If a fixed bounding box area is being used, then threshold the image plane within which the spot is located using a default sensitivity value (ex. 0.35). If doing this doesn't include the spot in the foreground, then increase the sensitivity value and repeat until at least the middle of the spot is included in the foreground. Following this, compute the centroid of the portion of the spot that is in the foreground. Finally, crop the original unthresholded image plane down to a box of the desired size (ex. 7x7 or 9x9); this will be the region of interest (ROI) for that spot.

3b. If a variable bounding box area is being used, then threshold the image plane within which the spot is located using a very low sensitivity value (0.02). Compute the area of the bounding box around the portion of the spot that is included within the foreground. If this area is outside of the desired range for the bounding box area, modify the sensitivity value and re-threshold. Stop either once a bounding box of the desired area is attained, or a certain number of iterations are performed. In the latter case, where no bounding box of appropriate area can be found, do not include this spot in any further calculations. Finally, crop the original unthresholded image plane down to a box equal to the bounding box, with an extra two (this number is controlled by the ebbr variable) pixels in every direction; this will be the region of interest (ROI) for that spot. 

4. Reshape the intensity values in the ROI into a single vector, for the purpose of the curve fitting.

5. Fit the intensity values in the ROI to the 2D Gaussian surface I = B+A*e^-((x-x0)^2/(2*sigmax^2) + (y-y0)^2/(2*sigmay^2)), where I is the intensity, B is the background value, A is the peak amplitude of the surface, x0 and y0 are the coordinates of the surface's peak, and sigmax and sigmay are the standard deviations in the x and y directions of the surface.

6. Extract the portion of the Gaussian that was a part of the original bounding box (i.e. if a variable bounding box size was used, exclude the extra two pixels that were added on each side prior to the curve fitting). Over this portion of the Gaussian, sum the intensity values and subtract the product of the background intensity and the number of pixels within this area (i.e. compute the total intensity above the background). 

7. Discard spots that either have a corresponding background intensity of less than zero (this happens occasionally if the spot that was proposed isn't actually isolated in the image), or that have an eccentricity of above 0.75. Eccentricity is calculated as sqrt(1 - min(sigmax^2 / sigmay^2 , sigmay^2 / sigmax^2)), and a very large eccentricity indicates an elliptical spot/ROI that is likely the result of multiple different NPCs being included within one ROI.

8. Display a histogram showing integrated intensity vs. the count of isolated spots. 

## III. Using the Program

a. Setting up the Configuration File

By default, this script relies upon a configuration file named SNSDconfig.txt, which is stored in the MATLAB home folder. A sample layout for this file can be found included in this repository.

The first line of the configuration file must be of the format imageNames = name1,name2,name3 etc, where name1 is the name (without file extension) of the first stack of images being fed into the program, and so on. The name to the left of the equals sign is not important, but the use of commas as delimiters and the lack of spaces after any of the names is. 

The second line should be of the format outputFileName = name, where name is the name of the file to which all data will be output at the conclusion of the script. Again, the spacing is important.

The third line should be of the format showFigures = x, where x is either 'none', 'some', or 'very' (without quotes). The value of x will determine the value of the verbose variable in the script. If 'none' is chosen, only the final summary histogram will be output by the script. If 'some' is chosen, then for each spot a heatmap of the Gaussian fit and a heatmap of the ROI will be displayed as well. If 'very' is chosen, then a number of other images (such as the binary mask of the entire uncropped image and of the ROI) will be displayed for each spot, which can be useful when troubleshooting.

The fourth line should be of the format boundingBoxAreaRange = x,y where x and y are both integers. If x is equal to y, then the script will use a fixed bounding box size that is equal to x. If this is desired, then x and y should be equal to an odd perfect square (either 49 or 81 seems to work best). If x is unequal to y, then the script will use a variable bounding box size, where x is the lowest size that the script will accept and y is the largest acceptable size. When using a variable box size, 16,36 seemed to work well - note that because two pixels will be added to each side of the box before fitting, the box will end up being between 8x8 and 10x10. 

Each subsequent line will correspond to the location, within one of the image stacks specified in the first line of the config file, of a single isolated NPC spot. Each of these lines will have the format p,x,y,i where p is the plane of the image in which the spot is located, x and y are the coordinates of the brightest pixel in the proposed spot within plane p, and i is the number of the image stack within which the spot is located. To be clear, if a spot is located in the 10th plane of the 3rd image stack that is listed in the first line of the config file, then the line corresponding to that spot should be 10,x,y,3. 

In future, there are a number of additional parameters that are currently hardcoded that should be moved to the configuration file, including the excess bounding box radius (ebbr) used when working with variable bounding boxes, and the maximum Gaussian eccentricity that is allowed for a given spot to be counted (currently set to 0.75). As a further note on the ebbr variable, always make sure that it is set to 0 before running the script using fixed bounding boxes (or just insert a line of code that sets ebbr to 0 automatically if fixed bounding boxes are used). 

b. Optimizing Parameters

The main parameters which could be tweaked to try and optimize results are the upper and lower limits for acceptable bounding box size, the excess bounding box radius that is employed when using variable bounding boxes, and the maximum Gaussian eccentricity that is allowed for a given spot to be counted. Another parameter (currently hardcoded) that could be modified is the default sensitivity value used to threshold images when in fixed bounding box mode, although the current value works well in every circumstance that has been tested. 

I have found that the most consistent results are obtained when the range of acceptable bounding box areas is 16 <= area <= 36 and an excess bounding box radius of 2 is used. Generally, these parameters should hold values that ensure that a) every noticably bright pixel of a spot is included within the initial unexpanded bounding box (as it is only this area that will later be integrated over) and b) the expansion of the bounding box by ebbr pixels in every direction adds pixels whose intensities are very close to the background while not adding pixels from other nearby NPC spots. This will allow the entirety of the Gaussian peak to be included within the final bounding box, while also providing a sample of the background that should allow it to be estimated more accurately by the fitting function. For more details, please refer to the Discussion of the Interim Project Report (Section IV. c. of this README). 

When it comes to the maximum permissible eccentricity and the default sensitivity value used to threshold images in fixed bounding box mode, there is no particular rationale underlying my choices of their current values. The former should only be modified if the script's output suggests that there are many single isolated spots being excluded by the eccentricity filter or many double/triple spots that are able to slip past it. The latter should only be modified if microscope images are taken in such different conditions from what I used that the default value in the script is too high and results in the entire image being included within the foreground.

c. Understanding the Output

The output of the script has two main components: a .csv file containing information about the Gaussian fit for each set of spot coordinates provided in the input. The column labels should largely be self-explanatory. Peak is the peak intensity at the centre of the Gaussian, while int is the total integrated intensity. The value in the usable column is 1 if a bounding box of appropriate area was found for the spot and it passed the eccentricity filter, and 0 otherwise. Only spots with a usable value of 1 are included in the summary histogram. Finally, the bounding box area column displays the final bounding box area that was calculated for that spot, including the pixels added by the excess bounding box radius expansion. This area may be extremely large in cases where the script was not able to find an appropriate bounding box.

The histogram is very standard. In theory, the distribution of spot intensities should itself be a Gaussian, centred on the 'true' fluorescence intensity of a single isolated NPC. In reality, there will be some error in the process by which the intensity is calculated, and there will also be some non-isolated spots that contain two NPCs that are not eliminated by the filtering process. Because of this, fitting the distribution of spot intensities with a Gaussian Mixture Model (k = 2) is recommended, as this can separate the single NPC spots from the double NPC spots. So long as the peaks of the resulting Gaussians are in approximately a 2:1 ratio, it can reasonably be concluded that this separation was successful.

d. Applications of the Output

When the script is run on many sets of images, each with a different NPC tagged, a standard curve of mNG intensity vs. stoichiometry can be constructed in order to estimate the intensity of a single mNG molecule. The resulting value should then be used as input to the NuclearIntensityIntegrator script, in order to calibrate its results by converting total nuclear intensity to a number of mNG molecules.
