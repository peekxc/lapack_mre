project('hello', 'cpp')
cpp = meson.get_compiler('cpp')
link_args = []

py_mod = import('python')
py = py_mod.find_installation(pure: true)
py_dep = py.dependency()
message('Python path =' + py.full_path())
message('Numpy version =' + run_command(py, ['-c', 'import numpy; print(numpy.__version__)'], check: true).stdout().strip())

## Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

## Recommended on Numpy docs: https://numpy.org/doc/stable/f2py/buildtools/meson.html
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()
inc_np = include_directories(incdir_numpy, incdir_f2py)

## TODO: fix this
use_scipy_wheels = true
if use_scipy_wheels
  ## Link statically to SciPy's openblas wheel, if available
  blas = dependency('scipy-openblas', required: true, method : 'pkg-config')
  inc_blas = include_directories(blas.get_variable(pkgconfig: 'includedir'))
  link_args += blas.get_variable(pkgconfig: 'extralib').split()
else
  blas = dependency('openblas', required: true, method : 'pkg-config', static: false)
  inc_blas = include_directories(blas.get_variable(pkgconfig: 'includedir'))
  link_args += blas.get_variable(pkgconfig: 'extralib').split()
  # inc_blas = include_directories('/usr' / 'local' / 'opt' / 'openblas' / 'include')
  # link_args += '-L/usr/local/opt/openblas/lib' # brew recommended links to compile
  # link_args += '-L/usr/local/lib' # Link to default system library on UNIX systems
endif

## It just doesn't make sense: I get library not loaded /usr/local/lib/libopenblas, but there is no explicit call to 
## link openblas via -lopenblas or anything; only -lopenblas_python, but that is provided by the scipy wheel. 
executable(
  'hello_lapack', 
  'hello.cpp', 
  link_args: link_args,
  dependencies: [blas], 
  include_directories: [inc_np, inc_blas] # ['/usr/local/lib/', '/usr/local/gfortran/lib/']
)