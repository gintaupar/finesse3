.. include:: /defs.hrst

.. _cameras_usage:

Capturing beam images with camera objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :mod:`finesse.detectors.camera` module includes a number of Camera objects for
probing single pixels, slices or full images of beam(s). These cameras are split into
categories according to these dimensions as well as another two categories corresponding
to whether it is detecting the field amplitude and phase (i.e. a sub-class of
:class:`.ComplexCamera`) or the pixel intensity (ie. a sub-class of
:class:`.CCDCamera`).

Equations describing the quantities computed by each of these two latter categories can
be found in :ref:`camera_equations`.

.. rubric:: Complex (field) cameras

Listed below are the camera objects which are used for probing the per-pixel field
amplitude and phase for a beam at a specific frequency.

.. autosummary::

    finesse.detectors.camera.FieldCamera
    finesse.detectors.camera.FieldScanLine
    finesse.detectors.camera.FieldPixel

These are listed in order of decreasing dimensions.

.. rubric:: CCD cameras

Listed below are the camera objects which are used for probing the per-pixel beam
intensity of all (non-audio) fields present at the detection node.

.. autosummary::

    finesse.detectors.camera.CCD
    finesse.detectors.camera.CCDScanLine
    finesse.detectors.camera.CCDPixel

These are listed in order of decreasing dimensions.

.. note::

    All camera classes take the optional argument ``w0_scaled`` in their constructors.
    This flag is True by default and indicates that all :math:`x` and :math:`y`
    coordinates should be scaled by the waist size of the beam. Any outputs (including
    plots) will then be in units of the beam waist. By setting this flag to False, you
    can prevent this behaviour and instead have the :math:`x` and :math:`y` axes in
    units of metres.

Detecting a single pixel
````````````````````````

A single pixel of a beam profile can be detected with :class:`.Pixel` type cameras. The
return type of this detector will be a single number (complex for a field, real for a
CCD). A simple example below shows how one can use a :class:`.FieldPixel` to detect, for
example, the central pixel (at :math:`x = y = 0`) of a beam distribution on reflection
from a mirror whilst tilting the mirror in the x-direction.

.. jupyter-execute::

    import finesse
    finesse.configure(plotting=True)

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1 L=1
    m ITM R=0.99 T=0.01
    gauss gL0 L0.p1.o q=(-2+3j)

    # a FieldPixel detector - defaults to x = y = 0 and f = 0
    fpx refl_px ITM.p1.o

    modes(maxtem=2)
    """
    )

    out = model.run("xaxis(ITM.xbeta, lin, 0, 100u, 100)")
    out.plot();

Note that, as implied in the code above, the constructor for :class:`.FieldPixel` uses
default values of zero for x and y (so that the central pixel is the default), and it
probes at a frequency offset of 0 Hz (carrier) by default too.

Probing a slice of a beam
-------------------------

Slices of a beam profile are detected with :class:`.ScanLine` type cameras. These
detectors return a 1D :class:`numpy.ndarray` where the scanning axis is specified at
construction time via the ``direction`` argument. An example is shown below, where the
x-axis of a beam image, in terms of the intensity, is computed and plotted using a
:class:`.CCDScanLine` detector.

.. jupyter-execute::

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1 L=1
    m ITM R=0.99 T=0.01
    gauss gL0 L0.p1.o q=(-2+3j)

    # A CCDScanLine detecting x in [-3, 3] with y = 0 by default
    ccdline refl_x ITM.p1.o xlim=3 npts=100

    modes(maxtem=2)
    """
    )

    out = model.run("noxaxis()")
    out.plot();

A single value can be specified for x, y when either of these is the scanning axis -
this will then create the scanning axis array ranging from, e.g., :math:`-|x|` to
:math:`+|x|` if :math:`x` is the axis to be scanned. The value for the non-scanned axis
will take on a default of zero if not specified.

Capturing a full beam image
```````````````````````````

Full two-dimensional beam profiles can be captured with :class:`.Image` type cameras.
Detectors of this type return a 2D :class:`numpy.ndarray`. An example is shown below,
using the same model as above but now taking the full beam distribution image in terms
of the intensity with a :class:`.CCD` detector.

.. jupyter-execute::

    # Remove the CCDScanLine detector from above
    model.remove(model.refl_x)

    # Add a CCD with x in [-1, 1] and y in [-1, 1] with 300 pixels per axis
    model.parse("""
    ccd refl_im ITM.p1.o xlim=1 ylim=1 npts=80
    """)

    out = model.run("noxaxis()")
    out.plot();


Producing animated beam images
``````````````````````````````

.. warning::

    Depending on your working environment, animations may not be enabled for displaying.
    If you are working in a Jupyter Notebook you may need to use `%matplotlib notebook`
    or `%matplotlib qt` to set up your environment, other wise a single fixed image
    will be displayed.

Animated plots of beam distributions can be created by adding a full image detector
(i.e. a :class:`.CCD` or :class:`.FieldCamera`) to a model, then scan over a model
parameter (using the ``xaxis`` command for example) to produce 3D data which will be
animated using the parameter scan axis when plotting via :meth:`.ArraySolution.plot`.

An example of this is shown below, where a distorted beam is injected into a Fabry-Perot
cavity. A cavity scan is then performed whilst a CCD captures the beam image, at each
mirror detuning, to produce a mode scan animation.

.. jupyter-execute::
    :hide-output:

    model = finesse.Model()
    model.parse(
    """
    l L0 P=1

    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-10
    s CAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=10

    cav FP ITM.p2.o

    tem(L0, n=0, m=1, factor=0.9)
    tem(L0, m=1, n=0, factor=0.9)
    tem(L0, n=0, m=2, factor=0.8)
    tem(L0, n=1, m=1, factor=0.4)
    tem(L0, n=2, m=0, factor=0.5)

    ccd trns ETM.p2.o xlim=2.5 ylim=2.5 npts=80

    modes(maxtem=2)
    """
    )

    out = model.run("xaxis(ITM.phi, lin, 0, 180, 60)")
    figures, animations = out.plot(cmap="hot")

.. jupyter-execute::
    :hide-code:

    from IPython.display import HTML
    HTML(animations["trns"].to_jshtml())
