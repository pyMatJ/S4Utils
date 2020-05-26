..  _Installation:

Download and installation
==========================

Here are the guidelines to install S\ :sup:`4` and the Python API either as a :ref:`Conda package <S4Conda>` or :ref:`Building from source <S4FromSource>`

------------------------------------------
Installing S\ :sup:`4` and the Python API
------------------------------------------

..  _S4Conda:

Conda package
-------------

If you run Anaconda, S\ :sup:`4` can be installed as a conda package. From an anaconda prompt

.. code-block:: sh
    
    conda config --add channels conda-forge  
    conda upgrade conda
    
To add the channel conda-forge as one of the default channel for downloading packages and upgrade all you packages. Create a new dedicated conda environment (*e.g* with python 3.7:

.. code-block:: sh

    conda create -n s4 python=3.7 
    
Activate the environment and install the package:

.. code-block:: sh

    conda activate s4  
    conda install -c paulgoulain s4  

You should be able to import the S\ :sup:`4` module in a python prompt::

    import S4 as S4

..  _S4FromSource:

Building from source
---------------------

You can download `S4 <https://github.com/victorliu/S4>`_  on github, however we recommend using slightly more up-to-date forks of the original code, by cloning *e.g.*

.. code-block:: sh
    
    git clone https://github.com/phoebe-p/S4 
    cd S4


or 

.. code-block:: sh

    git clone https://github.com/marcus-o/S4
    cd S4


Windows
~~~~~~~

To be added

Linux
~~~~~ 
Two options are proposed, depending on which optional libraries you use for the linear algebra (BLAS and Lapack), Fourier transform (FFT) and matrix factorization (CHOLMOD). 

The first one makes use of the (free but proprietary) `Math Kernel Library (MKL) from Intel <https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html>`_, which largely improves performances (at least on Intel CPUs). :ref:`jump <S4MKL>`

The second one uses only open-source libraries (`OpenBlas <http://www.openblas.net/>`_, `FFTW3 <http://www.fftw.org/>`_ and `CHOLMOD (SuiteSparse) <https://github.com/DrTimothyAldenDavis/SuiteSparse/releases>`_, at the risk of lower performances. :ref:`jump <S4Open>`

..  _S4MKL:
 
**Proprietary option**


The first step is to `install the MKL libraries <https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html>`_:

.. code-block:: sh

    sudo bash
    # <type your user password when prompted.  this will put you in a root shell>
    # cd to /tmp where this shell has write permission
    cd /tmp
    # now get the key:
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    # now install that key
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    # now remove the public key file exit the root shell
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    exit
    # add the repository
    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    # update the package list
    sudo apt update

You can now install one of the MKL releases. Since you have to specify the component you wish to install, you can `check for the version list <https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html>`_

.. code-block:: sh

    ## exemple intel-mkl-2020.0-088
    sudo apt install intel-mkl-<VERSION>.<UPDATE>-<BUILD_NUMBER> 


We also need the `OpenMP` and `suitespares-dev` libraries:

.. code-block:: sh
    
    # OpenMP
    sudo apt install libomp-dev
    # CHOLMOD
    sudo apt-get install libsuitesparse-dev 

The correct compilation flags have to be set in the `Makefile` to use the MKL libraries. Edit the following lines if necessary, the other options should work:

.. code-block:: sh

    BLAS_LIB = -lmkl_rt
    LAPACK_LIB = -lmkl_rt
    
    FFTW3_INC =
    FFTW3_LIB = -lfftw3

Then from a command prompt in the `S4` folder:

.. code-block:: sh
    
    make boost
    make S4_pyext
    
.. note:: You might have to adapt the end of the Makefile to use the correct ``pip`` installer. ``pip3`` is used by default, but if you want to install S\ :sup:`4` in a particular environment, you should first activate it before running `make S4_pyext`. 

In case an error occurs, it is most probably because the `MKL` libraries were installed somewhere where the compiler does not find them (*e.g* in ``/opt/``). They should be linked as (might require `sudo`):

.. code-block:: sh

    ln -s $(MKLdir)/*.so /usr/lib/$(where_to_store)
    
*i.e* in my case:

.. code-block:: sh

    sudo ln -s /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/*.so /usr/lib/x86_64-linux-gnu/

Everything should compile without any errors or warnings. 

..  _S4Open:

**Open source option** 


We start by installing the OpenBlas and CHOLMOD libraries:

.. code-block:: sh

    ## linear alegra blas and lapack
    sudo apt-get install libopenblas-dev 
    ## fourier transform fftw3 
    sudo apt-get install libfftw3-dev 
    ## CHOLMOD for some fourier decompositions
    sudo apt-get install libsuitesparse-dev 

Check that the Makefile contains the following lines (all other options should already be set):

.. code-block:: sh

    BLAS_LIB = -lblas
    LAPACK_LIB = -llapack
    
    FFTW3_INC =
    FFTW3_LIB = -lfftw3

Then compile:

.. code-block:: sh
    
    make boost
    make S4_pyext

.. note:: You might have to adapt the end of the Makefile to use the correct ``pip`` installer. ``pip3`` is used by default, but if you want to install S\ :sup:`4` in a particular environment, you should first activate it before running `make S4_pyext`. 

Everything should compile without any errors or warnings. 

------------------
Installing S4Utils
------------------

Installing S4Utils only requires ``pip``. From a command prompt in the main S4Utils directory:

.. code-block:: sh

    pip install .
    
Should run without any trouble. The last tutorial script :ref:`Tuto4-MIM` (`download script <../../examples/MIM_DispersiveGrating.py>`_) uses some of S4Utils functions and should thus run if the installation is successful. 
