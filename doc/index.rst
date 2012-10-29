.. batse5bp documentation master file, created by
   sphinx-quickstart on Thu Oct 25 17:09:59 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


..  |nbsp|  unicode:: 0xA0 
   :trim:

.. |---| unicode:: U+02014 .. em dash
   :trim:

The batse5bp package
====================

.. toctree::
   :maxdepth: 2

Introduction
------------

The batse5bp package provides persistant local access to BATSE GRB data hosted
at the CGRO SSC.  It accesses both the BATSE GRB **catalog** (a collection of
basic properties of all detected bursts, such as directions, intensities, and
durations), and **detailed GRB data** (e.g., time
series and spectral data).

The last official BATSE GRB catalog was the 4B catalog.  The "current" catalog
hosted at the SSC includes many GRBs observed subsequent to the 4B catalog.  A
5B catalog, including data spanning to the end of the CGRO mission, was
planned but never completed (in particular, the sky exposure and efficiency
tables were never computed).  The "5bp" part of the package name denotes "5B,
preliminary" (5Bp).

The full 5Bp GRB catalog is fetched, parsed, and stored locally when the
user first loads the catalog.  Subsequent loads access the stored copy.
Detailed GRB data, which would be voluminous if downloaded en masse, is
fetched transparently only when the user first attempts to access it.
Only the requested data element, for a particular GRB, is fetched; it is
cached locally for subsequent access.

Cached data is stored in a local file hierarchy with a structure resembling
that of the SSC archive.  Files and folders may be deleted within the local
filesystem when no longer needed; they will be re-fetched as necessary.  There
is no monolithic database file, and no dependencies on database software.


Installation
------------

``batse5bp`` is currently a pure-Python package with only ``numpy`` and
``pyfits`` as dependencies (some example code here also relies on
``matplotlib``).  It is provided as a gzipped tar file containing the source
code, "batse5bp-0.1.tar.gz".  It could be installed by unpacking and using the
standard ``python setup.py install`` command.  But the package is evolving so
the recommended method is to use ``pip``, which will allow easy uninstallation
and replacement:

.. code-block:: sh

    pip install batse5bp-0.1.tar.gz

To uninstall it:

.. code-block:: sh

    pip uninstall batse5bp

..
    If you do not have access to your computer's global ``site-packages``
    directory, or if you'd like to keep your ``batse5bp`` activity
    insulated from your other Python activity, use ``virtualenv`` to
    create a virtual Python environment, and then install within that
    environment with ``pip``.  If you haven't yet installed and set up
    ``virtualenv``, I recommend using ``virtualenv-burrito``, which
    will install both ``virtualenv`` and ``virtualenv-wrapper``, and
    also modify your shell startup script to properly set them up
    for use.

    * virtualenv-burrito_

    * `Getting Started with virtualenv and virtualenvwrapper in Python <http://www.silverwareconsulting.com/index.cfm/2012/7/24/Getting-Started-with-virtualenv-and-virtualenvwrapper-in-Python>`_ [ignore the *Installation* and *Configuration* sections; ``virtualenv-burrito``
    handles that]

    * `Quick & Dirty Virtualenv (& Virtualenvwrapper)`_

    .. _virtualenv-burrito:  https://github.com/brainsik/virtualenv-burrito


Loading the catalog
-------------------

To establish access to the 5Bp data, use the ``load_catalog()`` function::

    from batse5bp import load_catalog

    # Load the catalog using the default database location (in the CWD):
    grbcat = load_catalog()

The first time it's run, ``load_catalog()`` creates a database directory for
storing fetched data, grabs the 5Bp catalog files from the SSC, parses
them, collects the data in an object, and saves a serialized version of
the object in the database for future access.  The return value is the
object providing access to the data, described below.

By default, the database directory is named "BATSE_5Bp_Data" and is created or
accessed in the CWD.  To create it elsewhere, or with a different name, or
to access a previously loaded database stored elsewhere, provide
the database directory path as an argument to ``load_catalog()``.

The GRBCollection object
------------------------

``load_catalog()`` returns a ``GRBCollection`` object that stores the
data in ``GRB`` objects, one for each GRB.

A ``GRBCollection`` is a *mapping* (like a Python ``dict``), providing access to
the data for individual GRBs by using the BATSE trigger number (an integer) as
a key.  For example, ``grbcat[105]`` returns a ``GRB`` instance providing
access to the data for BATSE trigger 105, the first GRB observed by BATSE.

A ``GRBCollection`` is an *ordered* mapping, so if you iterate over it, you
will access the GRB instances in order of the BATSE trigger number.  (A
standard Python ``dict`` returns its keys and items in an arbitrary order.)

For convenience, besides the standard ``dict`` keyword
access, you can access GRB data via *attributes*, in two ways:

* ``grbcat.tN`` returns the GRB instance for trigger number *N*, e.g.,
  ``grbcat.105`` returns the instance for trigger 105.

* ``grbcat.bYYMMDD`` returns a **list** of GRB instances for bursts
  that occured on the specified date, in the YYMMDD format that is
  traditional for identifying GRBs.  For example, ``grbcat.b910421``
  returns a list with one element, the GRB instance with the data
  for trigger 105, which happens to be GRB |nbsp| 910421.  A list is
  returned because more than one burst
  may have occured on a given day.  (Although conventional GRB
  nomenclature distinguishes such bursts with lower case letter
  suffixes, the current BATSE catalog does not provide such full
  IDs, which would require coordination with catalogs from other
  missions.)

So all three of the following ``GRB`` instances are the same::

    grb = grbcat[105]
    grb2 = grbcat.t105
    grb3 = grbcat.b910421[0]


The GRB object: Catalog data
----------------------------

You will access the data and metadata for a GRB via the attributes
and methods of its ``GRB`` instance.

Suppose you load the catalog in an interactive IPython session, and
grab the ``GRB`` instance for trigger 105, as follows:

.. sourcecode:: ipython

    In [1]: from matplotlib.pylab import *

    In [2]: ion()  # for interactive plotting in examples below

    In [3]: from batse5bp import load_catalog

    In [4]: grbcat = load_catalog()
    Loaded summary data for 2702 GRBs comprising the 5Bp catalog.

    In [5]: grb = grbcat.t105

The string representation of ``grb`` (what gets shown by a ``print``
statement) provides a handy summary of the data in the burst's
catalog entries:

.. sourcecode:: ipython

    In [6]: print grb
    Trigger:  105
    Name:  4B_910421
    TJD:   8367.384766
    (RA, Dec) (long, lat) (deg):    (270.68   24.76)  ( 21.19   50.75)
    Drxn err (deg):    0.53
    Complete?  True
    Ch 1-4 fluence (ergs):  8.7e-07  1.3e-06  2.0e-06  1.0e-06
    Ch 1-4 S/N:             8.0e+01  9.2e+01  1.1e+01  5.8e+00
    64ms F_pk (cts), S/N, t:   12.76   23.46     3.84
    T50, T90:    1.79    5.18
    A:  Ulysses, PVO rate increase
    O:  Background source (Cyg X-1)
    Q:  Data gap during burst accumulation; HV turnoff at T+38 s

Each item can be accessed as an attribute of ``grb``; for example:

.. sourcecode:: ipython

    In [11]: grb.trigger
    Out[11]: 105

    In [12]: grb.desig
    Out[12]: '910421'

    In [13]: (grb.F64ms, grb.F64ms_err)
    Out[13]: (12.761, 0.544)

These attributes contain, respectively, the trigger number, the designator in
YYMMDD format (but omitting a possible lower-case suffix), and the 64 |nbsp|
ms time scale peak flux and its error.  See the ``grb`` module API
documentation (in the Module Index, linked below) for more information about
the attributes holding the catalog data.

Of note is the inclusion of information from the BATSE "Comments" table
in the ``GRB`` instance:

.. sourcecode:: ipython

    In [14]: grb.comments
    Out[14]: 
    [('A', 'Ulysses, PVO rate increase'),
     ('O', 'Background source (Cyg X-1)'),
     ('Q', 'Data gap during burst accumulation; HV turnoff at T+38 s')]

There were unusual circumstances for a number of GRBs that can impact
interpretation of the data; the comments describe most such cases.
A burst may have multiple comments associate with it, each with a
flag identifying the comment type.  The ``comments`` attribute is a
list of 2-tuples containing the flag and its associated comment.

To understand the comments and the other catalog attributes, the
authoritative reference is the catalog description text available at
the SSC.  The package provides local copies of the text for the
various tables comprising the catalog in the ``docn`` module.  For
example, to see an explanation of the comments:

.. sourcecode:: ipython

    In [15]: from batse5bp import docn

    In [16]: print docn.comments
    The COMMENTS Table contains comments relevant to BATSE gamma-ray burst data
    found in the GROSSC BATSE burst catalog. Each category of comment is sorted by
    a flag identification and ordered by BATSE trigger number. A gamma-ray burst
    may have more than one entry or may have no entry.

    Flag  Definition
    Q     comments on data quality
    A     additional observations by other instruments
    O     general comments
    L     comments on the gamma-ray burst coordinates
    T     comments on the gamma-ray burst duration

The ``docn`` module contains similar text strings for the three other tables
comprising the catalog:  ``basic`` (with the trigger number, data, direction,
overwrite flag, etc.), ``flux`` (with peak fluxes and fluences), and ``durn``
(with durations).

The ``docn`` module also provides four convenience functions that will
open important web pages at the SSC in your default browser; for
convenience, these are all inserted into the ``batse5bp`` namespace:

``docn.show_ssc()``
  Show the COSSC web site.

``docn.show_4B()``
  Show the 4B catalog web site.
  (The 4B catalog remains useful as the largest catalog with 
  exposure and efficiency tables; also, the 4B bursts have
  official designators, including letter suffixes in cases where
  multiple bursts happened on the same day.)

``docn.show_5Bp()``
  Show the 5Bp ("current") catalog web site.

``docn.show_problems()``
  Show the BATSE data problems archive.

The problems archive has text files identified by trigger number that
describe problems with the data for a few dozen bursts, beyond what
is mentioned in the Comments table.


The GRB object: Detailed data
-----------------------------

``GRB`` instances provide access to some of the detailed data for each
GRB stored at the SSC in FITS files and files of other formats.

Light curve plots
~~~~~~~~~~~~~~~~~

The simplest detailed data are light curve plots, stored as GIF image files at
the SSC.  Two plots are archived: one shows 64 |nbsp| ms count data from each
of BATSE's four discriminator (DISC) channels, the other shows the summed
light curve.  The ``show_lc()`` method fetches the images, creates a
single image from them, and displays it using your default image viewer.
The images are stored in the local database for future use.  This is a
good way to get a quick sense of what the data look like, but note that
the archived light curves often show only a subset of the data.

Compiled 64 ms count data (ASCII64)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The BATSE team has compiled three different DISC data types into a four-
channel "64 |nbsp| ms" binned counts data set for each burst, stored as an
ASCII file.  These data are summed count *rates* (float-valued) from the subset
of BATSE's eight detectors triggered by a burst.  Two of the data types
contain genuine 64 |nbsp| ms data associated with a trigger: the PREB (pre-
burst) and DISCSC (discriminator science) data.  To extend coverage further,
both before and after the time covered by PREB and DISCSC, the BATSE team used
the continuously-recorded 256 |nbsp| ms DISCLA data to estimate 64 |nbsp| ms
count rates prior to the PREB and after the DISCSC data.  The early ASCII data
thus has contiguous sets of 16 bins with identical count rates. Since the
trigger need not be at the boundary of a 256 |nbsp| ms DISCLA interval, there
is often a single set of fewer than 16 ASCII bins with identical rates, just
before the PREB data.

A detailed description of the ASCII 64 |nbsp| ms data, copied from the SSC,
is available as a string in the ``ascii64`` module:  ``ascii64.docn``.

*Note that the BATSE DISC pulse height channels are conventionally numbered
from 1 to 4, but are accessed in Python code using 0-based indexing, i.e.,
numbered from 0 to 3.*

``GRB`` instances have an ``ascii64`` attribute that provides access to the
ASCII 64 |nbsp| ms data via an ``ASCII64`` object with several attributes and
one method providing different views of the data.  The aim of this object is
to report integer-valued *counts* (not rates), in the actual 64 and 256 |nbsp|
ms bins used for measurements, by re-binning the subdivided DISCLA data,
and converting the rate estimates back to counts.
Presently this is only correctly done for the early DISCLA data.  Fortunately,
most bursts are short enough that late-time DISCLA data will typically be
ignored.  But for long bursts, the user should use this ASCII data with caution,
finding a way to handle the late-time subdivided DISCLA data.  For the early
mission, the DISCSC data duration is about 4 |nbsp| m; later in the mission
the flight software was reconfigured to provide ~10 |nbsp| m of DISCSC data.
The late-time DISCLA bins appear after these intervals.

For the attribute descriptions below:

* Counts are in integer arrays, and times are in float arrays.

* Unless otherwise specified, the time arrays list bin *boundaries*,
  consecutively, so if there are N bins, the time array will have length N+1.

* Times are in seconds relative to the trigger time.

* Count data is provided in 4-row arrays, so that ``counts[i,j]`` is the
  number of counts in time bin ``j`` for channel ``i``, and ``counts[2,:]``
  is an N-vector containing the count time series for channel 2.

One set of attributes makes the data with different binning separately
available.  ``early`` and ``t_early`` provide the counts and
time bin boundaries for the early 256 |nbsp| ms DISCLA data. 
Similarly, ``late`` and ``t_late``
provide counts and times for the 64 |nbsp| ms PREB and DISCSC
data (but including any late-time subdivided DISCLA data).
In cases where there is a single truncated DISCLA bin before the
PREB data, the ``mid`` and ``t_mid`` attributes contain the single
photon count value in the truncated bin (as a length-1 NumPy integer array),
and its time boundaries (as a length-2 float array).  If there
is no truncated DISCLA bin, the value of ``mid`` is ``None``.

This block of code prints the last few counts in the ``early`` array, the
``mid`` counts (if present, else ``None``), and the first few counts in the
``late`` array::

    asc = grb.ascii64
    print 'Early:', asc.early[0,-3:]
    print 'Mid:', asc.mid
    print 'Late: ', asc.late[0,:3]


A second set of attributes concatenates the data into single count and time
arrays.  ``counts`` contains all of the counts in time sequence.  ``times``
contains the boundaries of all of the time bins.

The ``rates`` method returns a 3-tuple of NumPy arrays built from the
concatenated data::

    (centers, widths, rates) = grb.ascii64.rates()

``centers`` contains the center times of the bins, ``widths`` contains
the durations of the bins, and ``rates`` is a 4-row array containing counting rate
estimates found by dividing the counts by the bin width.  This is
intended for plotting light curves.  The following block of code
uses ``matplotlib`` to plot the light curves for each channel in a common figure::

    asc = grb.ascii64
    t, w, r = asc.rates()
    figure(figsize=(12, 7))
    for ch in range(4):
        plot(t, r[ch,:], '-', lw=2, alpha=.6)
    axvline(c='k')  # show the trigger time
    xlabel(r'$t$ (s from trigger)')
    ylabel(r'$Rate$ (counts/s)')

(This plots the light curve as a curve; a histogram would be a more
accurate representation of the data.)

Detector response matrices (DRMs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``GRB`` instances have a ``discsc_drms`` attribute providing access to the
DRMs for the DISC channels of triggered detectors, via a ``DRMs_DISCSC``
object.  This object contains ``DRM`` objects that hold the DRM data for each
triggered detector.  It also provides access to a DRM that sums the response
of the triggered detectors.

Users should note that the archived DRMs are fairly coarsely sampled.  Also, the
response for DISC channel 4 is relatively flat over the tabulated range and is
simply truncated at the high-energy end of the range (about 2 GeV).  Integrals
of the product of an incident flux model and the channel response will be
inaccurate if the incident flux does not fall with energy at least as quickly as
1/E**2.

``discsc_drms.n_det`` is the number of triggered detectors.  

Two ``discsc_drms`` attributes describe the common format
of the DRMs for the triggered detectors and the summed response:

``n_ch``
  The number of pulse height channels (always 4 for DISC data).  

``n_E``
  Each channel's response is described by a vector evaluating the response
  function at ``n_E`` values of incident photon energy (at least in *theory*
  |---| see the ``E_bins`` attribute, below).

A DRM is a matrix defined so ``drm[i,j]`` is the response of channel ``i`` for
the ``j``'th value of incident photon energy.

``discsc_drms.detectors`` is a list of ``DRM`` instances, one for each
triggered detector.

``DRM`` instances have the following attributes:

``det_num``
  The detector number, an integer from 1 to 8.

``ch_bins``
  The pulse height bin boundaries, in units of nominal deposited energy
  in keV.

``E_bins``
  The ``n_E+1`` incident photon energy boundaries ("edges" in BATSE parlance) defining
  the DRM.  Peculiarly, this is provided as if some kind of binning was done in
  incident photon energy, making the energy assignment for a particular DRM
  entry ambiguous.  Currently I suggest using the center of the ``E_bins`` bin; this
  needs further exploration.

``drm``
    The DRM as an ``n_ch`` by ``n_E`` array.

``crv``
    A list of ``n_ch`` channel response vectors, each an array of
    length ``n_E``.  These are just different views of the DRM,
    useful for calculating expected count rates for individual channels.

``start``
    The response is vanishingly small at low energies.  ``start[i]``
    is the index of the first value of ``j`` for which ``drm[i,j]`` is
    non-zero.


The *summed* response of the triggered detectors is available via
similar attributes of ``discsc_drms`` itself:
``ch_bins``, ``E_bins``, ``drm``, ``crv``, and ``start``.

The following block of code plots the summed response for all four DISC
channels in a single figure.  Note that ``E_vals`` is set equal to the
incident energy bin centers.

::

    drms = grb.discsc_drms
    d0 = drms.detectors[0]
    d1 = drms.detectors[1]

    if 1:
        E_vals = (drms.E_bins[:-1] + drms.E_bins[1:]) / 2.
        for i in range(drms.n_ch):
            semilogx(E_vals, drms.crv[i], '-', lw=2)
        xlabel(r'$E$ (keV)')
        ylabel(r'$R_i(E)$ (cm$^2$)')

The ``DRMs_DISCSC`` provides a capability to numerically integrate the product
of an incident spectrum model and the response function for each channel; the
result is the expected number of detected photons ("counts") in the channel.
The calculation is done using an interpolatory inner product quadrature rule
that uses evaluations of the response function and the spectrum at *different*
points, to account for different availability and scales of variation of the
factors.  Currently, only one rule is available:  For each incident energy
interval in a channel response vector, the rule uses the two response values at
the interval boundaries, and three spectrum values within the interval (located
at the zeroes of a 2nd-degree Legendre polynomial over the interval).

Current quadrature support is only for the *summed* response.

To set up the rules for the summed response channels for a burst's DISC data,
call the ``set_sum_prodquad`` method of the burst's ``DRMs_DISCSC`` object.
This defines a product quadrature rule for each channel; it also compiles
the nodes (incident energy values used in the rules) and makes them
available as an array via the ``quad_nodes`` attribute, so
``quad_nodes[i,j]`` returns the ``j'th`` energy value for channel ``i``.
The ``chan_quad(i, spec)`` method calculates the quadrature for channel
``i`` for input spectrum ``spec``.  The spectrum can be provided as a
function (of incident photon energy), or as a vector of values pre-calculated
on the nodes (this will be faster if results are sought for more than one
channel).  If a spectrum function is used that requires arguments besides
the photon energy, those arguments can be provided to ``chan_quad`` as
additional arguments; e.g., for a spectrum with parameters ``p1`` and
``p2``, the signature would be: ``chan_quad(i, spec, p1, p2)``.

The following block of code pre-calculates spectrum values and calculates
the expected number of counts for BATSE channels 1 and 2 (indices 0 and 1).
It assumes a spectrum function ``spec(E)`` exists, and that it can accept
array arguments.

::

    drms = grb.discsc_drms
    drms.set_sum_prodquad()
    svals = spec(drms.quad_nodes[0,:])  # channels have common E values
    exp_cts_0 = drms.chan_quad(0, svals)
    exp_cts_1 = drms.chan_quad(1, svals)


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

