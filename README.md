# MASNUM-ocean-circulation-model
The case about the climatological simulation for the Indian Ocean.

The MASNUM ocean circulation model is a two-time-level, three-dimensional numerical ocean circulation model with a two-level, single-step Eulerian forward-backward time-differencing scheme. A mathematical model of large-scale oceanic motions was based on the terrain-following coordinated, Boussinesq, Reynolds-averaged primitive equations of ocean dynamics. A simple but very practical Eulerian forward-backward method was adopted to replace the most preferred leapfrog scheme as the time-differencing method for both barotropic and baroclinic modes. The forward-backward method is of second-order of accuracy, computationally efficient by requiring only one function evaluation per time step, and free of the computational mode inherent in the three-level schemes.

makefile: Compiling the code.
masnum3d.f90: The main code.
netcdf_mod.f90: Outputting and inputting the netcdf file in convenience.
time_mod.f90: Converting the time variables.
