import os 
import importlib

def scipy_build_config():
  import scipy as sp
  config = sp.show_config(mode="dicts")
  m_info = config['Machine Information']
  b_info = config['Build Dependencies']
  return b_info

if __name__ == '__main__':

  print("Detecting scipy openblas")
  ob32_found = importlib.util.find_spec("scipy_openblas32") is not None
  ob64_found = importlib.util.find_spec("scipy_openblas64") is not None
  if not ob32_found and not ob64_found:
    raise ModuleNotFoundError("Did not find scipy openblas wheels. Please 'pip install scipy_openblas{32|64}' first")
  blas_variant = '64' if ob64_found else '32'

  ## For now, we just write to the current working directory
  basedir = os.getcwd()
  openblas_dir = os.path.join(basedir, ".openblas")
  pkg_config_fname = os.path.join(openblas_dir, "scipy-openblas.pc")
  os.makedirs(openblas_dir, exist_ok=True)

  ## Get SciPy build configuration
  b_info = scipy_build_config()

  ## Write the pkgconfig file
  openblas = importlib.import_module(f"scipy_openblas{blas_variant}")
  pkg_config = openblas.get_pkg_config().split('\n')
  for i, config_str in enumerate(pkg_config):
    if config_str[:5].upper() == 'LIBS:' and len(config_str) <= 6:
      config_str += r"-L${libdir}"
      pkg_config[i] = config_str
    elif config_str[:8].upper() == 'EXTRALIB':
      extra_lib = pkg_config[i]
      pkg_config[i] += f" -L{b_info['blas']['lib directory']}"
      pkg_config[i] += f" -L{b_info['lapack']['lib directory']}"
    elif config_str[:6].upper() == 'CFLAGS':
      pkg_config[i] += f" -I{b_info['blas']['include directory']}"
      pkg_config[i] += f" -I{b_info['lapack']['include directory']}"

  pkg_config = '\n'.join(pkg_config)
  with open(pkg_config_fname, "wt", encoding="utf8") as fid:
    fid.write(pkg_config.replace("\\", "/"))
  
  ## Print status and finish
  print(f"Finished writing '{pkg_config_fname}' to directory '{openblas_dir}' ")