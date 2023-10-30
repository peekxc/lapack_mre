To build, activate a `virtualenv` and then use: 

> python gen_openblas_pc.py

to generate a pkgconfig file suitable for usage with meson. By default, a `pc` file is created in generated `.openblas` directory yielding the linkage information for `scipy-openblas{32|64}`, depending on which are installed. 

From the project directory, use the following command to create a meson build directory: 

> PKG_CONFIG_PATH=.openblas meson setup --wipe build 

And then compile everything with the following:

> cd build && meson compile 

The final executable tests a few basic functionalities of BLAS/LAPACK. 

> ./hello_lapack