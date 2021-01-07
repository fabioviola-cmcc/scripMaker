# scripMaker.py

## Basic usage

Invoke the script with:

```
$ python scripMaker.py --grid=grid.nc --mask=mask.nc --output=out.nc
```

## Requirements

scripMaker.py can be run in an environment configured like this:

```
$ conda create -n scrip
$ conda activate scrip
$ conda install "python<=3.9.0a0"
$ conda install netCDF4
```

The full list of packages at the end of the installl process should be like:

```
$ conda list
# packages in environment at /users_home/opa/resm-dev/.conda/envs/scrip:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
blas                      1.0                         mkl  
bzip2                     1.0.8                h7b6447c_0  
ca-certificates           2020.12.8            h06a4308_0  
certifi                   2020.12.5        py38h06a4308_0  
cftime                    1.3.0            py38h6323ea4_0  
curl                      7.71.1               hbc83047_1  
hdf4                      4.2.13               h3ca952b_2  
hdf5                      1.10.4               hb1b8bf9_0  
intel-openmp              2020.2                      254  
jpeg                      9b                   h024ee3a_2  
krb5                      1.18.2               h173b8e3_0  
ld_impl_linux-64          2.33.1               h53a641e_7  
libcurl                   7.71.1               h20c2e04_1  
libedit                   3.1.20191231         h14c3975_1  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 9.1.0                hdf63c60_0  
libgfortran-ng            7.3.0                hdf63c60_0  
libnetcdf                 4.7.3                hb80b6cc_0  
libssh2                   1.9.0                h1ba5d50_1  
libstdcxx-ng              9.1.0                hdf63c60_0  
mkl                       2020.2                      256  
mkl-service               2.3.0            py38he904b0f_0  
mkl_fft                   1.2.0            py38h23d657b_0  
mkl_random                1.1.1            py38h0573a6f_0  
ncurses                   6.2                  he6710b0_1  
netcdf4                   1.5.3            py38hbf33ddf_0  
numpy                     1.19.2           py38h54aff64_0  
numpy-base                1.19.2           py38hfa32c7d_0  
openssl                   1.1.1i               h27cfd23_0  
pip                       20.3.3           py38h06a4308_0  
python                    3.8.5                h7579374_1  
readline                  8.0                  h7b6447c_0  
setuptools                51.0.0           py38h06a4308_2  
six                       1.15.0           py38h06a4308_0  
sqlite                    3.33.0               h62c20be_0  
tk                        8.6.10               hbc83047_0  
wheel                     0.36.2             pyhd3eb1b0_0  
xz                        5.2.5                h7b6447c_0  
zlib                      1.2.11               h7b6447c_3  
```

## Error codes

The following error codes are raised in case of exceptions:

1. One or more input parameters are missing
2. Grid file not found
3. Mask file not found
4. Opening of the output file failed

More will be added soon.