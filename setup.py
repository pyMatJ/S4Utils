from setuptools import setup, find_packages

setup(
   name='S4Utils',
   version='0.0.1',
   description='Utilities function for the Python API of S4',
   author='Mathieu Jeannin, Paul Goulain',
   author_email='',
   package_dir={"": "src"},
   packages=[''],
   install_requires=['numpy', 'scipy', 'matplotlib', 'pyvista'] #external packages as dependencies
)