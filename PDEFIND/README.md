This folder contains the code to perform the PDE learning aspect of our study.

Run **script_PDEFIND_properror_sf_pruning.py** to run the PDEFIND algorithm implemented in the study.

**script_PDEFIND_properror_sf_pruning.py**--scripting version of PDE_find_properror_sf_pruning.ipynb to be able to run on servers

**PDE_find_properror_sf_pruning.ipynb** - ipython notebook to run PDE-FIND implementations. In the second cell, the variable "comp_str"specifies which denoising strategy one wishes to use (nn, finite_differences, splines (meaning cubic bisplines), NCV_bisplines (meaning local cubic bisplines with a GLS error model), or global_NCV_bisplines_3 (meaning global cubic bisplines with a GLS error model) ) and the variable "model_str" specifies which model one wants to consider (diffadv, fisher, fisher_nonlin). The third cell then updates other various aspects of the study, as detailed throughout the paper.

**Properror analyze results.ipynb** - Once one has performed the PDE-FIND calculations and saved results in the folder "pickle_data", they can plot their results by using this ipython notebook. In the second cell here, the list "model_str_list" specifies which mathematical models will be considered. 

**nonlin_fisher_IP.ipynb** - ipython notebook to perform the inverse problem methodology discussed in Section 3(e) of Lagergren et al (2020).

**make_learned_eqn_tables.py**- Create the tables of learned equations provided in our supplementary material

**PDE_FIND2.py** -- code used to implement the PDE-FIND algorithm

**pickle_data/**-- contains output from the PDEFIND algorithm

**Data/**--data produced from the data denoising part of the study

**Learned equations tables/**--learned equation results

