from setuptools import setup, find_packages

setup(
   name='S4Utils',
   version='0.0.1',
   description='Utilities function for the Python API of S4',
   author='Mathieu Jeannin, Paul Goulain',
   author_email='',
   package_dir={"S4Utils": "S4Utils"},
   packages=['S4Utils'],
   install_requires=['numpy', 'scipy', 'matplotlib'], ### pyvista/vtk dependency not supported yet
   ##, 'pyvista'], #external packages as dependencies
   package_data={'S4Utils': ['S4Utils/examples/*']},
)
