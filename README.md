# Supplementary Material for manuscript "Ray tracing-based optimization of a 45 kW solar-thermal Cassegrain dish"
## By Ramón Pujol-Nadal, Luis Guerreiro, and Gabriel Cardona

Contents of the repository:

* `1.Cassegrain_designs`: Folder containing the Cassegrain designs.
  * `create_scene_cassegrian.py`: Script for generating the geometries of the Cassegrain concentrator, based on six geometric parameters, for its analysis in OTSunWebApp.
  * `optic_materials.zip`: Optical materials of the constructed designs (in OTSunWebApp file format) for the simulations.
 
* `2.Optical_efficiency`: Folder containing files related to the optical efficiency described in Section "3.2. Optical efficiency".
  * `results_set_1.txt.gz`: Output of the simulation of the first set of designs (compressed).
  * `results_set_2.txt.gz`: Output of the simulation of the second set of designs (compressed).
  * `plots_results_set_1.py`: Script used to generate Fig. 9, which shows optical efficiencies for the first set of designs. Generates the file `Optical_Efficiency_set_1.png`.
  * plots_results_set_2.py: Script used to generate Fig. 10, which shows optical efficiencies for the second set of designs. Generates the file `Optical_Efficiency_set_2.png`.
    
* `3.Power_balance`: Folder containing files related to the power balance analysis described in Section "3.3. Power balance".
  * `Power_balance.py`: Script used to generate the values shown in Table 5, which lists the ten designs with the highest absorbed power by the receiver. Generates the file `Power_balance.png` (not shown in the paper).

* `4.Tracking_error`: Folder containing files related to the tracking error described in Section "3.4. Tracking error".
  * `results_effiopt_vs_angle_T.txt`: Output of the simulation of Cases 1 and 2 in the transversal plane.
  * `results_effiopt_vs_angle_L.txt`: Output of the simulation of Cases 1 and 2 in the longitudinal plane.
  * `plots_tracking_error.py`: Script used to generate Fig. 11, which shows optical efficiency in transversal and longitudinal planes. Generates the file `tracking_error.pdf`.

* `5.Flux_receiver`: Folder containing files related to the flux radiation on the receiver described in Section "3.5. Radiation distribution on the receiver".
  * `simulation_outputs`: Folder with the output files from the simulations of Cases 1 and 2 using the Spectral Analysis mode in OTSunWebApp (some files therein are compressed).
  * `ASTMG173-direct.txt`: Spectral direct solar radiation.
  * `flux_distribution.py`: Script used to...
  * `plots_tracking_error.py`: Script used to generate Figs. 12-13, which show radiation flux distribution for Cases 1 and 2. Generates the figures `heat_map_3D_Case_1.pdf`, `heat_map_3D_Case_2.pdf`, `heat_map_Case_1.pdf`, and `heat_map_Case_2.pdf`.
    
* `6.Appendix`: Folder containing files related to the Appendix.
  * `simulation_outputs`: Folder with output files for each number of emitted rays: 10k, 50k, 100k, 200k, 500k, and 1M. 
  * `N_ray_error_sigmas.py` and `N_ray_error_sigmas_b.py`: Scripts used to obtain Tables A.1 and A.2. Also generate the figures `plot_errors.pdf` and `plot_errors_sigma.pdf`.
