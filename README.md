# U-Pb-LA-ICP-MS-Reduction
## Program for reducing U-Pb data collected by laser ablation. Sister to Lasertram-DB (Lubbers et al., 202x)

## Main To-Do's:

1. Add File upload button

	* Hopefully do without secondary trigger
	* Need to get everything else talking to the uploaded file. Currently, some functions depend on data import via pandas
	* Setup callback in class with Watch=True that send upload to input_data widget possibly
    
2. Add errors for all ratios, not just 206/238

	* Use ratio_buttons widget as input into get_regressions, get_residuals, get_approved_regressions
	* Will require some string formatting
	* Need to expand list in get_approved function as well.
  
