# An_2021_DBS_TRD_sim

* These codes are about building individual virtual brains (surface-based models) and simulating response networks evoked by local deep brain stimulation, reflecting fiber tract engagement. For more details on this, see An et al., 2021, Neuroimage.

* The codes are composed of a python-based main script that performs network simulation and matlab-based scripts that derive the inputs necessary to conduct the simulation, and the execution order is as follows. It is assumed that structural brain reconstruction using T1-weighted, diffusion-weighted MRI and CT scans of individual subjects has been carried out beforehand (refer to the tvb reconstruction pipeline for this task).



  * 1 `surf_const.m`    

    This code excludes the brain stem and very few vertices located at unknown positions from the reconstructed cortical/subcortical structures.



  * 2 `conn_gen.m`    

    This code generates connectivity (connection strength matrix and track-length matrix) for the surface-based model composed of ~10000 nodes. This work is performed by connecting the nodes at both ends based on the geometry of each of the 100 000 fiber tracts regardless of brain parcellation.



  * 3 `volt_dist_cal.m`    

    This code calculates the voltage distribution around the stimulation position based on the finite difference method.



  * 4 `stim_gen.m`    

    This code derives the stimulus input matrix to be applied to each node based on the voltage distribution calculated according to the stimulation position, and the geometries of the fiber tracts. This work is performed by applying a stimulus to the nodes at both ends of the activated fiber tracts, assuming that the fiber tract is activated when the voltage across it exceeds a certain threshold.



  * 5 `surface_based_sim.ipynb`    

    This code simulates the brain response characteristics induced by local stimulation, i.e., by selective stimulation of fiber tracts, by applying the inputs obtained above.
