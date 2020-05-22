Tutorials
=========

Tutorial 1: A single, flat interface
-------------------------------------

In this example, we compute the reflection and transmission coefficients at a single, flat interface between a supestrate (incident medium) and a substrate. The simulation can be run either in two or three dimensions, to demonstrate the syntax adjustments in creating the S\ :sup:`4` simulation object.


Once S4 has been installed, you should be able to import it, along with numpy and matplotlib.pyplot::

    import S4 as S4
    import numpy as np
    import matplotlib.pyplot as plt

As stated in the `S4 documentation <https://web.stanford.edu/group/fan/S4/units.html>`_, "S\ :sup:`4` solves the linear Maxwell’s equations, which are scale-invariant. Therefore, S\ :sup:`4` uses normalized units, so the conversion between the numbers in S4 and physically relevant numbers can sometimes be confusing." In order to clearly link the reduced units and the actual physical ones, we keep track of the speed of light (assumed to be 1 in S\ :sup:`4`)::

    c_const = 3e8

and the general scale is defined by the period of the lattice::

    a = 1 ## period, normalized units. This defines the length scale

So if we (arbitrarily) decide that `a` here is in units of micrometers, we can run the simulation at a wavelength of 1 micron by setting the following::

    lbda = 1e-6 ## say we compute at 1µm wavelength
    f = c_const/lbda ## ferquency in SI units
    f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]

We pick the dimension of the simulation (2 or 3), then we set the lattice period and the number in-plane Fourier expansions order to be used (`S4.New <https://web.stanford.edu/group/fan/S4/python_api.html#S4.New>`_)::

    Dimension = 3 # 2 or 3

    if Dimension == 2:
        ### 2 dimensional simulation, 1D lattice
        S = S4.New(Lattice = a, 
                    NumBasis = 1)
    elif Dimension == 3:
        ### 3 dimensional simulation, 2D lattice
        S = S4.New(Lattice = ((a,0), (0,a)), 
                    NumBasis = 1)

In the case of a 2D simulation, the lattice is one-dimensional and has a simple period `a`. For the 3D simulation, the lattice dimension is two and we choose a square lattice 
of period `a` in both dimensions. Since the layer is unpatterned, we only need to compute the 0\ :sup:`th` plane wave order, and hence ``NumBasis=1`` is enough.

We then setup the main parameters of the simulation: the angle of incidence :math:`\theta`, and the refractive indices of the superstrate and substrate::

    theta = np.arange(0,89,2) ## angle of incidence in degree

    n_superstrate = 1 ## incident medium index
    eps_superstrate = n_superstrate**2 ## incicent medium permittivity
    n_substrate = 2 ## substrate index
    eps_substrate = n_substrate**2 ## substrate permittivity 

S\ :sup:`4` uses (possibly complex and anisotropic) permittivities instead of refractive index, which is used when declaring the materials used in the simulation (`S.SetMaterial <https://web.stanford.edu/group/fan/S4/python_api.html#S4.Simulation.SetMaterial>`_)::

    S.SetMaterial(Name='Air', Epsilon = eps_superstrate)
    S.SetMaterial(Name='Substrate', Epsilon = eps_substrate)

The names used in declaring the materials (``'Air'`` and ``'Substrate'``) will be used to set the permittivities of the layers that we declare, setting their names and thickness (`S.AddLayer <https://web.stanford.edu/group/fan/S4/python_api.html#S4.Simulation.AddLayer>`_)::

    AirThick = 1 ## in reduced units, same scale as a
    SubThick = 1 ## in reduced units, same scale as a

    S.AddLayer(Name='Air', Thickness=AirThick, Material='Air')
    S.AddLayer(Name='Substrate', Thickness=SubThick, Material='Substrate')

The layers are added in order, from the top (incident layer) to the bottom (substrate) one. 

At this point the simulation is set up, and we can now run it to retrieve the reflection and transmission coefficients. We initialize the arrays storing the results::

    Rte = np.zeros(len(theta))
    Rtm = np.zeros(len(theta))
    Tte = np.zeros(len(theta))
    Ttm = np.zeros(len(theta))

and set the (single) frequency at which the simulation runs (`S.SetFrequency <https://web.stanford.edu/group/fan/S4/python_api.html#S4.Simulation.SetFrequency>`_)::

    S.SetFrequency(f0)

The bulk of the calculation is in the main loop over the angle of incidence::

    for ii, thi in enumerate(theta):
        ########
        ## TM excitation
        S.SetExcitationPlanewave(
                IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
                sAmplitude=0.,
                pAmplitude=1.,
                Order=0)
        inc, r = S.GetPowerFlux(Layer='Air',zOffset=0)
        fw, _ = S.GetPowerFlux(Layer='Substrate', zOffset=SubThick)
        Rtm[ii] = np.abs(-r/inc)
        Ttm[ii] = np.abs(fw/inc)
        
        #########
        ## TE excitation
        S.SetExcitationPlanewave(
                IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
                sAmplitude=1.,
                pAmplitude=0.,
                Order=0)
        inc, r = S.GetPowerFlux(Layer='Air',zOffset=0)
        fw, _ = S.GetPowerFlux(Layer='Substrate', zOffset=SubThick)
        Rte[ii] = np.abs(-r/inc)
        Tte[ii] = np.abs(fw/inc)

The exciting wave properties are set through the `S4.SetExcitationPlanewave <https://web.stanford.edu/group/fan/S4/python_api.html#S4.Simulation.SetExcitationPlanewave>`_ method. The arguments are the couple of incidence angles (theta,phi) in spherical coordinates which represent the polar angle and the azimuthal angle. (Note that the S4 documentation uses the opposite denomination), the amplitude of the `s`- and `p`-components of the electric field, and the order (defaults to 0, see doc). Note that here we use the ``enumerate`` function from python which allows us to pass directly the local variable ``thi`` to the plane wave function. 

We use the method `S.GetPowerFlux <https://web.stanford.edu/group/fan/S4/python_api.html#S4.Simulation.GetPowerFlux>`_ which returns the integral of the `z`component 
of the Poynting vector over unit cell surface normal to the `z` direction, decomposed into forward and backward propagating modes. Hence, the incident and reflected power is obtained using::

    inc, r = S.GetPowerFlux(Layer='Air',zOffset=0)

so the calculation is lead in the first layer. The ``zOffset`` parameter specifies the vertical offset (in reduced units) from the *front* surface of the layer. This only matters for lossy layers. Conversly, the transmitted power is obtained with::

    fw, _ = S.GetPowerFlux(Layer='Substrate', zOffset=SubThick)

where we discard the value of the backward propagating wave from the substrate, since it should be 0 by definition. Note that we specified a ``zOffset`` equal to the layer thickness only for the sake of clarity. 

As stated in the `S4 documentation <https://web.stanford.edu/group/fan/S4/units.html>`_, the incident power is only unity at normal incidence. Hence, the values extracted from the simulation have to be normalized by the input power:: 

    Rtm[ii] = np.abs(-r/inc)
    Ttm[ii] = np.abs(fw/inc)

(note the sign convention for the reflected power) and we take the absolute value ``abs`` to cast the (otherwise complex) value returned by the simulation to a real number. 

We can plot the results and compare to the Fresnel coefficients formula, as well as the Brewster angle:

.. image :: Fresnel_plot.png

where the symbols are the results of the simulations, and the solid lines the analytical results.


Tutorial 2: A single slab, or Fabry-Pérot etalon
------------------------------------------------

We now quickly move on to a similar, simple example demonstrating the calculation of the transmission and reflection from a single dielectric slab, either as a function of angle of incidence or as a function of frequency. 
The begining of the script is exactly the same, setting the base unit and the dimensionality of the simulation::

    a = 1 ## period, normalized units. This defines the length scale
    Dimension = 3 # or 3
    if Dimension == 2:
        ### 2 dimensional simulation, 1D lattice
        S = S4.New(Lattice = a, 
                    NumBasis = 1)
    elif Dimension == 3:
        ### 3 dimensional simulation, 2D lattice
        S = S4.New(Lattice = ((a,0), (0,a)), 
                    NumBasis = 1)

and then the permittivities, material objects, thickness and layer objects::

    n_superstrate = 1 ## incident medium index
    eps_superstrate = n_superstrate**2 ## incicent medium permittivity

    n_slab = 2 ## substrate index
    eps_slab = n_slab**2 ## substrate permittivity 

    S.SetMaterial(Name='Air', Epsilon = eps_superstrate)
    S.SetMaterial(Name='Slab', Epsilon=eps_slab)

    AirThick = 1
    SlabThick = 1

    S.AddLayer(Name='AirTop', Thickness=AirThick, Material='Air')
    S.AddLayer(Name='Slab', Thickness=SlabThick, Material='Slab')
    S.AddLayer(Name='AirBottom', Thickness=AirThick, Material='Air')

Now two sets of calculations are possible. In the first one, we sweep over the angle of incidence as in the previous example. To do so, we first define a frequency at which we wish to run the calculation. Again, assuming that ``a`` above is in micrometers, we can define `e.g` a 1 micron wavelength and the corresponding frequencies, both in SI units (`f`) and reduced units(`f0`) and the angular range::

    lbda = 1e-6 ## say we compute at 1µm wavelength
    f = c_const/lbda ## ferquency in SI units
    f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]
    theta = np.arange(0,89,1)
    S.SetFrequency(f0)
    
The rest is exactly the same as above, and we obtain the angular reflectivity and transmittivity:

.. image :: FabryPerot_AnglePlot.png

In the second set of calculation, we fix the angle of incidence and set up a spectral range over which we run the calculation::

    theta = 15 ## fixed angle of incidence

    lbda = np.linspace(400e-9,1e-6,200) ## 400nm to 1µm
    f = c_const/lbda ## ferquency in SI units
    f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]

While the bulk of the code is the same, note that now we must remember to include the ``S.SetFrequency`` method in the main loop::
    
    for ii, fi in enumerate(f0):
        S.SetFrequency(fi)
        ########
        ## TM excitation
        S.SetExcitationPlanewave(
                IncidenceAngles=(theta,0), #S4 names are reversed (phi,theta)
                sAmplitude=0.,
                pAmplitude=1.,
                Order=0)
        inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
        fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
        Rtm[ii] = np.abs(-r/inc)
        Ttm[ii] = np.abs(fw/inc)
        
        #########
        ## TE excitation
        S.SetExcitationPlanewave(
                IncidenceAngles=(theta,0), #S4 names are reversed (phi,theta)
                sAmplitude=1.,
                pAmplitude=0.,
                Order=0)
        inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
        fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
        Rte[ii] = np.abs(-r/inc)
        Tte[ii] = np.abs(fw/inc)
    
which quickly allows to plot the results:

.. image :: FabryPerot_SpectrumPlot.png
