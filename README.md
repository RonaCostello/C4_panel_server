## A little tool for viewing relevant C4 expression for an orthogroup

### To use ###

#### First make sure the relevant packages are installed ####

There is a conda env file in this directory, use this to load up the environment:

  conda create --name your_env_name --file conda_env.txt
  conda activate your_env_name

If you use a different package manager, then just have a look at the top of the python script to see what you need. This should be okay, the packages are pretty standard ones.

#### Next identify the orthogroup of interest ####

For a given gene of interest, identify the orthogroup it belongs to using orthogroup_targetP_files

`grep -r 'AT5G19140' data/orthogroup_sequences/.`

#### Add the orthogroups you're interested in to a text file ####

New line for each orthogroup, no headings and put in the same folder as the python script

An example list is given - orthogroup_list_example.txt

#### Run the script! Add the orthogroup list file you made after calling the script ####

`python Panel_figure_C4_gene_expression.py your_list_of_orthogroups.txt`

#### Figures ####

The figures will save in the figures directory



If something gets changed and breaks, just let me know as this is a git repo so it can be recovered
