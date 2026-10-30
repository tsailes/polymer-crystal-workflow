# Build configuration only; official scientific sources are unmodified.
compiler = '/opt/pe_dft/envs/pe-gpaw-20261020/bin/mpicxx'
mpi = True
nolibxc = False
noblas = False
fftw = True
libvdwxc = True
scalapack = True
elpa = False
gpu = False
libraries = ['xc', 'blas', 'lapack', 'fftw3', 'scalapack', 'vdwxc']
include_dirs = ['/opt/pe_dft/envs/pe-gpaw-20261020/include']
library_dirs = ['/opt/pe_dft/envs/pe-gpaw-20261020/lib']
runtime_library_dirs = ['/opt/pe_dft/envs/pe-gpaw-20261020/lib']
extra_compile_args += ['-O3', '-fopenmp']
extra_link_args += ['-fopenmp']
