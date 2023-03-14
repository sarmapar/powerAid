<h2>Welcome!</h2>

Hi-C Poweraid is designed to help plan differential Hi-C experiments. <br> It is hosted at
[http://phanstiel-lab.med.unc.edu/poweraid](http://phanstiel-lab.med.unc.edu/poweraid).

This app consists of 2 main tabs, each with interactive plots that include hover effects for more detailed information and pan & zoom features. 

<h4>Power across Sequencing Depth</h4>
This tab is useful for determining the ideal total sequencing depth and number of replicates for an experiment and to investigate how differing dispersions can affect that decision. After choosing your desired values for dispersion, power threshold, and replicates, a plot will populate with the percentage of well-powered loops at various sequencing depths, where a well-powered loop is defined as one with a power at or above the selected power threshold to detect a 2-fold change in looping interactions. Multiple values for dispersion and replicates can be selected for easy comparison.

<h4>Power Across Loop Sizes</h4>
This tab provides more granular information on one set of parameters at a time. Here, you can make fine-tuned selections for fold change, dispersion, replicate number, and sequencing depth per replicate. The plots provided can then be used to determine power as a function of loop size. The top plot summarizes these into percentages per 10Kb binned distance, and the bottom plot shows the relative number of total and well-powered loops at each 10Kb binned distance. Here, a well-powered loop is defined as one with a power at or above the selected power threshold to detect the selected fold change in looping interactions. A summary of the results and chosen parameters is also displayed under the plots.


<h3> Running Hi-C Poweraid Locally </h3>
To run the app in your own R session, run the line

`runGitHub('poweraid', 'sarmapar', subdir = "poweraid_app")`

For this to work, make sure the package shiny is installed and loaded. All other packages will install and load on their own. It may take a moment the first time you run it, then it should be quick the subsequent times you launch it.

