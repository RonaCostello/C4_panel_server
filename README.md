## A little tool for viewing relevant C4 expression for an orthogroup

### To use ###

#### First make sure the relevant packages are installed ####

There is a conda env file in this directory, use this to load up the environment:

  `conda create --name your_env_name --file conda_env.txt
  conda activate your_env_name`

If you use a different package manager, then just have a look at the top of the python script to see what you need. This should be okay, the packages are pretty standard ones.

#### Next identify the orthogroup of interest ####

For a given gene of interest, identify the orthogroup it belongs to using orthogroup_targetP_files

`grep -r 'AT5G19140' data/orthogroup_sequences/.`

#### Add the orthogroups you're interested in to a text file ####

If there is only one orthogroup of interest, ignore this step:

Make a text file with each orthogroup on a new line. No headings and put it in the same folder as the python script

An example list is given - orthogroup_list_example.txt

#### Run the script! Add the orthogroup list file you made after calling the script ####

For a single orthogroup:

`python Panel_figure_C4_gene_expression.py OG00XXXXX` 

Or for multiple orthogroups:

`python Panel_figure_C4_gene_expression.py your_list_of_orthogroups.txt`

#### Optional: Specify a minTPM, this may be useful for very large orthogroups ####

Run the script, same as above but specify a minTPM for genes to be included in the plot and legend. The deafault is 0 (i.e. include all data points)

e.g. to set a minTPM of 100:

`python Panel_figure_C4_gene_expression.py OG00XXXXX 100`

#### Figures ####

The figures will save in the figures directory



If something gets changed and breaks, just let me know and I can recover a working version from github
