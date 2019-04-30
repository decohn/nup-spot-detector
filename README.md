# nup-spot-detector

Code to compute the integrated intensity of isolated Nup foci, designed for the Reyes-Lamothe Lab.

## Table of Contents

I. **Project Goals**

II. **Outline of Algorithm**

III. **Using the Program**

	a. Setting up the Configuration File (mention setting ebbr to 0 if using fixed bounding box size)
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

3b. If a variable bounding box area is being used, then threshold the image plane within which the spot is located using a very low sensitivity value (0.02). Compute the area of the bounding box around the portion of the spot that is included within the foreground. If this area is outside of the desired range for the bounding box area, modify the sensitivity value and re-threshold. Stop either once a bounding box of the desired area is attained, or a certain number of iterations are performed. In the latter case, where no bounding box of appropriate area can be found, do not include this spot in any further calculations. Finally, crop the original unthresholded image plane down to a box equal to the bounding box, with an extra two pixels in every direction; this will be the region of interest (ROI) for that spot. 

4. Reshape the intensity values in the ROI into a single vector, for the purpose of the curve fitting.

5. Fit the intensity values in the ROI to the 2D Gaussian surface I = B+A*e^-((x-x0)^2/(2*sigmax^2) + (y-y0)^2/(2*sigmay^2)), where I is the intensity, B is the background value, A is the peak amplitude of the surface, x0 and y0 are the coordinates of the surface's peak, and sigmax and sigmay are the standard deviations in the x and y directions of the surface.

6. Extract the portion of the Gaussian that was a part of the original bounding box (i.e. if a variable bounding box size was used, exclude the extra two pixels that were added on each side prior to the curve fitting). Over this portion of the Gaussian, sum the intensity values and subtract the product of the background intensity and the number of pixels within this area (i.e. compute the total intensity above the background). 

7. Discard spots that either have a corresponding background intensity of less than zero (this happens occasionally if the spot that was proposed isn't actually isolated in the image), or that have an eccentricity of above 0.75. Eccentricity is calculated as sqrt(1 - min(sigmax^2 / sigmay^2 , sigmay^2 / sigmax^2)), and a very large eccentricity indicates an elliptical spot/ROI that is likely the result of multiple different NPCs being included within one ROI.

8. Display a histogram showing integrated intensity vs. the count of isolated spots. 

## III. Using the Program
