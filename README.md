# Code developed and maintained exclusively by Eitan Rapaport, Andreas Benedikter, Michael Aivazis, Mark Simons 
# Planetary-Exploration-Tool (PET)

- Ground swaths calculated through look angle geometry instead of maximum incidence angle
- For ground intersects, nadir is sweeped instead of plane perpendicular to the velocity
- No SAR unwrapping error included
- Only works for repeat ground track orbits, otherwise we would have to do georeferencing to match swaths
- Plotting projections assuming biaxial ellipsoids
- Todo: Add creating 3D time-series from 1D time-series data

Tool for creating a digital twin of a mission to a solar system body.

# Installation Instructions

1. Have/install at least Python 3.7.2
2. Clone https://github.com/aivazis/mm.git
3. Create mm config file:

    a. Go to home directory

    b. Go to or create a .config directory

    c. Create a directory called "mm"

    d. Copy and paste the following into a file named "config.mm" in the "mm" directory,
**changing the python version to the one you are using**:
    
        # -*- Makefile -*-
        #

        # external dependencies
        # system tools
        sys.prefix := ${CONDA_PREFIX}

        # gsl
        gsl.dir := $(sys.prefix)

        # hdf5
        hdf5.dir := ${sys.prefix}
        hdf5.parallel := off

        # libpq
        libpq.dir := ${sys.prefix}

        # python
        python.version := 3.12
        python.dir := $(sys.prefix)

        # pybind11
        pybind11.dir = $(sys.prefix)

        # numpy
        numpy.dir := $(sys.prefix)/lib/python$(python.version)/site-packages/numpy/_core

        # pyre
        pyre.dir := $(sys.prefix)

        # install the python packages straight where they need to go
        builder.dest.pyc := $(sys.prefix)/lib/python$(python.version)/site-packages/

        # control over the build process
        # set the python compiler so we don't depend on the symbolic link, which may not even be there
        python3.driver := python$(python.version)
        
        # end of file

4. Create the mm yaml file:

    a. Go to home directory

    b. Go to .config directory

    c. Create a directory called "pyre"

    d. Copy and paste the following into a file named "mm.yaml" in the "pyre" directory, 
**changing the python version to the one you are using**::
    
        # -*- yaml -*-
        #

        # mm configuration
        mm:

          # targets
          target: opt, shared

          # compilers
          compilers: gcc, nvcc, python/python3

          # the following two settings get replaced with actual values by the notebook
          # the location of final products
          prefix: "{pyre.environ.CONDA_PREFIX}"
          # the location of the temporary intermediate build products
          bldroot: "{pyre.environ.HOME}/tmp/builds/mm/{pyre.environ.CONDA_DEFAULT_ENV}"
          # the installation location of the python packages
          pycPrefix: "lib/python3.12/site-packages"

          # misc
          # the name of GNU make
          make: make
          # local makefiles
          local: Make.mmm

        # end of file
        
5. Create a conda/mamba environment for the package

    a. Install conda/mamba if necessary

    b. Make a file in any directory and call it "pet.yaml"

    c. Copy and paste the following into the file:
    
        # -*- yaml -*-
        #


        name: PETenv

        channels:
          - conda-forge

            
        dependencies:
          - python
          - pyre
          - git
          - gcc
          - gxx
          - make
          - curl
          - gsl
          - hdf5
          - libpq
          - pip
          - setuptools
          - matplotlib
          - numpy
          - pip
          - pip:
            - cspyce
          - pybind11
          - pyyaml
          - ipython
          - pandas
          - matplotlib
          - tqdm
          - scipy
          - h5py
          - cartopy
          - xarray
          - alphashape
          - netcdf4

        # end of file

    d. Run the command "conda env create -f pet.yaml"
    e. Activate the environment

6. Go to the PET directory and run the command "python3 [PATH]/mm/mm.py" replacing [PATH] with the path to the mm directory 8
7. Ready to go!
