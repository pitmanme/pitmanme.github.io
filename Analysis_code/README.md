These are scripts 'academic code' and software for improved automation of some analysis 
and plotting tools. 

## Contained in this directory:

### For contact analysis:

* plot_contacts.py: driver code for calculating contacts for multiple simulations 
   type with multiple trials. 

* contact_analysis.py: functions for plot_contacts.py

### Running time series metric plotting:

* plot_rmsd_subplots.py: plots subplots for each type of simulation with multiple trial 
   data included. Includes a running average and raw data. 

* plot_avg_rmsd.py: plots average time series data over mutiple trials 
   (or from multiple simulations). 
   
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
