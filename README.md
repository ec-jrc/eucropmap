This work is in progress and may be subject to changes

# eucropmap

This repository contains the code developped for the paper:
"*From parcel to continental scale -- A first European crop type map based on Sentinel-1 and LUCAS Copernicus in-situ observations*" 
by Raphaël d'Andrimont, Astrid Verhegghen, Guido Lemoine, Pieter Kempeneers, Michele Meroni, Marijn van der Velde.

## Pre-processing
 `1_Legend_polygons_COPERNICUS_reorganise_v7` : Reclassify the Copernicus LUCAS polygons to match the EU crop map legend
 
 `2_Polygons_extract_GEE_reorganize_v7_S1`: Create a single .csv file with the S1 VV, VH and CR time series extracted for the LUCAS Copernicus polygons

  `AgregateS1timeSeries.R` : aggregate the S1 VV, VH from `S1_point_allV7_10days_10m_1Jan-31Dec_EU_ratio-db.csv` per group and compute statistic to prepare FIG 3 and FIG S7.
  
  
### Sentinel-1
 `GEE_S1_compositing.js` : creation of the 10-days VV and VH composite with Google Earth Engine
 `GEE_S1_polygons_extraction.js` : extraction of 10-days VV, VH and VH/VV time series for the Copernicus polygons (for feature selection)

### Classification
#### Parametrization
Parametrization of the four RF models 
#### Classification at scale
classification at scale

`eucropmap-classif-strata-v7.R`: The creation of a single geotiff from the different tiles and strata

### Post-processing
`maskeucropmap.py`: masking the EU crop map with the AW30D3, JRC GHSL European Settlement map 2015, JRC Global Water Surface product, Corine Land Cover 2018 and reprojecting to ETRS89-LAEA

### LUCAS

`CreateLUCAS2018validationPoints.R` : Creation of LUCAS validation point


## Accuracy 

`Validation_1_ClassifLucasPoints.R` : Validation from the LUCAS data (computing of FSCORE at the country level)

`Validation_2_GSAA.R` : Validation by comparing with the GSAA (i.e.farmers' declarations) 

`Validation_3_StatisticsNUTS2.R` : Validation by comparing Eurostat statistics



## Generate the figures of the paper

 `eucropmap_figures_tables_manuscript.R` to generate  FIG 3, 4, 5, 8 , 9 and  10 ; Supp Fig S1, S5, S6, S7
 

## Accuracy - pixac

Spatial accuracy assessment presented at IGARSS 2021

Verhegghen, A., D’Andrimont, R., Waldner, F., & Van der Velde, M. (2021). Accuracy assessment of the first EU-wide crop type map with LUCAS data. IEEE International Geoscience and Remote Sensing Symposium (IGARSS), 2021.
