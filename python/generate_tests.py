import sys, os
import f90nml
import subprocess

EQUILIBRIA = ['Adiabatic homogeneous',
              'Resistive homogeneous',
              'Gravito-acoustic waves',
              'Resistive tearing modes',
              'Resistive tearing modes with flow',
              'Suydam cluster modes',
              'Kelvin-Helmholtz',
              'Rotating plasma cylinder',
              'Kelvin-Helmholtz and current driven',
              'Rotating theta pinch']

NAMES = ['01-adiabatic-homo',
         '02-resistive-homo',
         '03-gravito-acoustic',
         '04-resistive-tearing',
         '05-resistive-tearing-flow',
         '06-suydam-cluster',
         '07-kh',
         '08-rotating-plasma-cyl',
         '09-kh-cd',
         '10-rotating-theta-pinch']

if sys.platform == 'darwin':
    os_suffix = '_MACOS'
else:
    os_suffix = '_LINUX'


if __name__ == '__main__':
    os.chdir("..")

    nml = {'gridlist': {}, 'equilibriumlist': {}, 'savelist': {}, 'filelist': {}}

    for idx, equil in enumerate(EQUILIBRIA):
        # gridlist
        nml['gridlist']['gridpoints'] = 101

        # equilibriumlist
        nml['equilibriumlist']['equilibrium_type'] = equil
        nml['equilibriumlist']['boundary_type'] = 'wall'

        # savelist
        nml['savelist']['write_matrices'] = False
        nml['savelist']['write_eigenvectors'] = False
        nml['savelist']['write_eigenfunctions'] = False
        nml['savelist']['show_results'] = False
        nml['savelist']['show_matrices'] = False
        nml['savelist']['show_eigenfunctions'] = False

        # filelist
        nml['filelist']['savename_eigenvalues'] = NAMES[idx] + os_suffix

        f90nml.write(nml, "parfile_test.par", force=True)

        subprocess.check_call(['./legolas', '-i', 'parfile_test.par'])

    os.remove("parfile_test.par")




