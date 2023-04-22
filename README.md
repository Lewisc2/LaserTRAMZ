# U-Pb-LA-ICP-MS-Reduction
## Program for reducing U-Pb data collected by laser ablation. Sister to Lasertram-DB (Lubbers et al., 202x)

## Main To-Do's:

1. Add File upload button

	* Hopefully do without secondary trigger
	* Need to get everything else talking to the uploaded file. Currently, some functions depend on data import via pandas
	* Setup callback in class with Watch=True that send upload to input_data widget possibly
	
	*Completed a work around that is much faster than the file upload button. Requires copy+pasting filepath*
    
2. Add errors for all ratios, not just 206/238

	* Use ratio_buttons widget as input into get_regressions, get_residuals, get_approved_regressions
	* Will require some string formatting
	* Need to expand list in get_approved function as well.
	
	*Completed 12/2022*
  
3. Add option to use 1st or 2nd order regression in second half of program titled 'CommonPbandnorm'
	*Completed*
	
4. Add Box and Whisker plot + option to plot up sample numbers for fliers / discordant data.
	* Points to be written down by user and inspected manually
	
	*Box and Whisker Plot Complete, but needs some reworking aesthetically*

5. In the very far future, try to calibrate 207/235 on the laser very well, then calculate 238 from 235 in order to deal with detector linearity issues for very high U zircons.
	* Cross check by calculating 207/235 on very old zircons - likely oracle, FC1, 91500, Tan-Bra
	* May need to add an option in second part of program to get ages from Wetherill Concordia. Tara-Wasserburg not so good for old ages.
