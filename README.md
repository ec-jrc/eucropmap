# eucropmap

This repository contains the code developped for the paper:
"*From parcel to continental scale -- A first European crop type map based on Sentinel-1 and LUCAS Copernicus in-situ observations*" 
by Raphaël d'Andrimont, Astrid Verhegghen, Guido Lemoine, Pieter Kempeneers, Michele Meroni, Marijn van der Velde.


## Pre-processing

### Sentinel-1

### LUCAS

`CreateLUCAS2018validationPoints.R` : Creation of LUCAS validation point


## Accuracy 

`ValidationClassifLucasPoints.R` Validation from the LUCAS data (computing of FSCORE at the country level)

Validation by comparing with farm declaration

Validation by comparing Eurostat statistics

## Generate the figures of the paper

 `eucropmap_figures_manuscript.R` to generate  FIG 3, 4, 5, 8 , 9 and  10 ; Supp Fig S1, S5, S6, S7
