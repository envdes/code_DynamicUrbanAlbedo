# code_DynamicUrbanAlbedo

[![DOI](https://zenodo.org/badge/770688781.svg)](https://zenodo.org/doi/10.5281/zenodo.10903399)

## Introduction

This repository is supplementary to the paper "**Improving Urban Climate Adaptation Modelling in the Community Earth System Model (CESM) Through Transient Urban Surface Albedo Representation**". 

The objectives of this project are:

- Modify CESM source code to realize transient urban albedo representation;
- Apply the new scheme for quantifying urban albedo cooling effects;
- Use simulation results for urban climate adaptation.

## Scripts and data

### [1_code_modification](./1_code_modification)

The standard source code comes from [CTSM](https://github.com/ESCOMP/CTSM), with the release tag: **[clm5.0.30](https://github.com/ESCOMP/CTSM/releases/tag/release-clm5.0.30)**.

- Create a new F90 file:
  
  - [src/biogeophys/UrbanDynAlbMod.F90](.//1_code_modification/src/biogeophys/UrbanDynAlbMod.F90)

- Apply the new module to:
  
  - [src/biogeophys/UrbanAlbedoMod.F90](./1_code_modification/src/biogeophys/UrbanAlbedoMod.F90)
  
  - [src/biogeophys/UrbanParamsType.F90](./1_code_modification/src/biogeophys/UrbanParamsType.F90)

- Add new name-list for case run:
  
  - [src/main/clm_varctl.F90](./1_code_modification/src/main/clm_varctl.F90)
  - [src/main/controlMod.F90](./1_code_modification/src/main/controlMod.F90)
  - [src/main/clm_driver.F90](./1_code_modification/src/main/clm_driver.F90)

- Add new name-list for case build:
  
  - [bld/CLMBuildNamelist.pm](./1_code_modification/src/bld/CLMBuildNamelist.pm)
  - [bld/namelist_files/namelist_defaults_clm4_5.xml](./1_code_modification/src/bld/namelist_defaults_clm4_5.xml)
  - [bld/namelist_files/namelist_definition_clm4_5.xml](./1_code_modification/src/bld/namelist_definition_clm4_5.xml)

- Add new name-list for initialization:
  
  - [src/main/clm_instMod.F90](./1_code_modification/src/main/clm_instMod.F90)

### [2_simulation_output_analysis]()

The scripts listed below are used for processing simulation output and visualization.

| Num. | Subject                                                      | Simulation                                                   | Output data process                                          | Visualization                                                |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 2.1  | [Roof albedo impacts on urban heat islands](./2_simulation_output_analysis/2.1_roof_albedo_impacts_UHI) | CNTL, ROOF_0.9, ROOF_DA                                      | Use [Export.ipynb](./2_simulation_output_analysis/2.1_roof_albedo_impacts_UHI/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.1_roof_albedo_impacts_UHI/Figure.ipynb) |
| 2.2  | [Roof albedo impacts on urban heat stress and indoor temperature](./2_simulation_output_analysis/2.2_roof_albedo_impacts_UHS) | CNTL, ROOF_0.9, ROOF_DA                                      | Use [Export.ipynb](./2_simulation_output_analysis/2.2_roof_albedo_impacts_UHS/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.2_roof_albedo_impacts_UHS/Figure.ipynb) |
| 2.3  | [Roof albedo impacts on surface energy budget](./2_simulation_output_analysis/2.3_roof_albedo_impacts_energy) | CNTL, ROOF_0.9, ROOF_DA                                      | Use [Export.ipynb](./2_simulation_output_analysis/2.3_roof_albedo_impacts_energy/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.3_roof_albedo_impacts_energy/Figure.ipynb) |
| 2.4  | [Urban surface heterogeneity in temperature](./2_simulation_output_analysis/2.4_urban_surface_heterogeneity_temperature) | CNTL, ROOF_DA, IMPROAD_DA, WALL_DA                           | Use [Export.ipynb](./2_simulation_output_analysis/2.4_urban_surface_heterogeneity_temperature/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.4_urban_surface_heterogeneity_temperature/Figure.ipynb) |
| 2.5  | [Urban landunit heterogeneity in temperature](./2_simulation_output_analysis/2.5_urban_landunit_heterogeneity_temperature) | CNTL, ROOF_DA, IMPROAD_DA, WALL_DA                           | Use [Export.ipynb](./2_simulation_output_analysis/2.5_urban_landunit_heterogeneity_temperature/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.5_urban_landunit_heterogeneity_temperature/Figure.ipynb) |
| 2.6  | [Spatial variation](./2_simulation_output_analysis/2.6_spatial_variation) | CNTL, ROOF_DA, IMPROAD_DA, WALL_DA                           | Use the 2040 outputs                                         | [Figure.ipynb](./2_simulation_output_analysis/2.6_spatial_variation/Figure.ipynb) |
| 2.7  | [Building energy in latitude](./2_simulation_output_analysis/2.7_building_energy_latitude) | CNTL, ROOF_DA, IMPROAD_DA, WALL_DA, ROOF_IMPROAD_DA, ROOF_IMPROAD_WALL_DA | Use the 2040 outputs                                         | [Figure.ipynb](./2_simulation_output_analysis/2.7_building_energy_latitude/Figure.ipynb) |
| 2.8  | [Building energy balance](./2_simulation_output_analysis/2.8_building_energy_balance) | CNTL, ROOF_DA, IMPROAD_DA, WALL_DA, ROOF_IMPROAD_DA, ROOF_IMPROAD_WALL_DA | Use [Export.ipynb](./2_simulation_output_analysis/2.8_building_energy_balance/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./2_simulation_output_analysis/2.8_building_energy_balance/Figure.ipynb) |

### [3_illustration](./3_illutration)

The figures listed below are used to illustrate details of the transient urban albedo scheme in CLMU.

| Subject                                                      | Visualization                                 |
| ------------------------------------------------------------ | --------------------------------------------- |
| Urban representation and parameterization in CLMU            | [Figure](./3_illustration/clmu.pdf)           |
| New module functionality                                     | [Figure](./3_illustration/dynalb.pdf)         |
| Comparing workflow between the default scheme and the new scheme | [Figure](./3_illustration/compare_scheme.pdf) |
| Reflectivity over urban surfaces                             | [Figure](./3_illustration/surface.pdf)        |
| Reflectivity over urban landunits                            | [Figure](./3_illustration/landunit.pdf)       |

### [4_supplimentary_information](./4_supplimentary_information)

The scripts listed below are used to show supplementary information such as input data, atmosphere variables, and computational performance.

| Num. | Subject                                                      | Analysis                                                     | Visualization                                                |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 4.1  | [Default input parameters in CLMU](./4_supplimentary_information/4.1_urban_parameter) | Use [Calculation.ipynb](./4_supplimentary_information/4.1_urban_parameter/Calculation.ipynb) to calculate urban fraction and roof albedo | [Figure.ipynb](./4_supplimentary_information/4.1_urban_parameter/Figure.ipynb) |
| 4.2  | [Actual surface albedo](./4_supplimentary_information/4.2_actual_surface_albedo) | Compare actual surface albedo                                | [Figure.ipynb](./4_supplimentary_information/4.2_actual_surface_albedo/Figure.ipynb) |
| 4.3  | [Regression](./4_supplimentary_information/4.3_regression)   | Use [Export.ipynb](./4_supplimentary_information/4.3_regression/Export.ipynb) to get the grid-cell level albedo and heat fluxes for [Calculation.ipynb](./4_supplimentary_information/4.3_regression/Calculation.ipynb) | NA                                                           |
| 4.4  | [Urban surface heterogeneity in energy budget](./4_supplimentary_information/4.4_urban_surface_heterogeneity_energy) | Use [Export.ipynb](./4_supplimentary_information/4.4_urban_surface_heterogeneity_energy/Export.ipynb) to get *.csv from 2015 to 2099 | [Figure.ipynb](./4_supplimentary_information/4.4_urban_surface_heterogeneity_energy/Figure.ipynb) |
| 4.5  | [Yearly atmosphere variables](./4_supplimentary_information/4.5_yearly_atmosphere_var) | Use [Export.ipynb](./4_supplimentary_information/4.5_yearly_atmosphere_var/Export.ipynb) to get the annual-mean outputs from 2015 to 2099 | [Figure.ipynb](./4_supplimentary_information/4.5_yearly_atmosphere_var/Figure.ipynb) |
| 4.6  | [Computational performance](./4_supplimentary_information/4.6_computational_performance) | Compare [timing log](./4_supplimentary_information/4.6_computational_performance/timing_log) | [Figure.ipynb](./4_supplimentary_information/4.6_computational_performance/Figure.ipynb) |

## Acknowledgments

- This work used the [ARCHER2 UK National Supercomputing Service](https://www.archer2.ac.uk). 
  The authors would like to acknowledge the assistance given by Research IT and the use of the HPC Pool and Computational Shared Facility at The University of Manchester. 
- The support of [Douglas Lowe](https://github.com/douglowe) and Christopher Grave from Research IT at The University of Manchester is gratefully acknowledged. 
- [Zhonghua Zheng](https://github.com/zhonghua-zheng) appreciates the support provided by the academic start-up funds from the Department of Earth and Environmental Sciences at The University of Manchester.
- [Yuan Sun](https://github.com/YuanSun-UoM) is supported by the PhD studentship of Zhonghua Zheng's academic start-up funds.
- Contributions from [Keith W Oleson](https://staff.ucar.edu/users/oleson) are based upon work supported by the NSF National Center for Atmospheric Research, which is a major facility sponsored by the U.S. National Science Foundation under Cooperative Agreement No. 1852977.
- Lei Zhao acknowledges the support of the U.S. National Science Foundation (CAREER award Grant 2145362).
- The authors declare no conflict of interest.
