This software improves automation of some analysis and plotting tools. This directory will grow with more tools as I add features. 

## Contained in this directory:

### For contact analysis:

* plot_contacts.py: driver code for calculating contacts for multiple simulations 
   type with multiple trials. 

* contact_analysis.py: functions for plot_contacts.py

[Example contact analysis plot output](https://github.com/pitmanme/pitmanme.github.io/blob/master/Analysis_code/example_outputs/avg_contacts.pdf)

### Running time series metric plotting:

* plot_rmsd_subplots.py: plots subplots for each type of simulation with multiple trial 
   data included. Includes a running average and raw data. 

* plot_avg_rmsd.py: plots average time series data over mutiple trials 
   (or from multiple simulations). 
   
[Example RMSD subplot plot output](https://github.com/pitmanme/pitmanme.github.io/blob/master/Analysis_code/example_outputs/peptide_rmsd_subplots.pdf)

[Example average RMSD analysis plot output](https://github.com/pitmanme/pitmanme.github.io/blob/master/Analysis_code/example_outputs/Peptide_avg_rmsd.pdf)
   
## Requirements:

* python 3.11.0
* pandas 1.5.3
* mdanalysis  2.4.2
* matplotlib 3.7.1

Other package versions may work but I have not verified them. 

## To run:

* Place code in the same directory (or in a searchable location on your machine). 
* For contact analysis, edit input variables in _plot_contacts.py_ and then run with:
```
          python plot_contacts.py
```
* For time series metric plotting (such as RMSD), edit input variables in 
  _plot_rmsd_subplots.py_ or _plot_avg_rmsd.py_ and run with
```
          python plot_rmsd_subplots.py
          python plot_avg_rmsd.py
```
