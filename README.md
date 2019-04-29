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


