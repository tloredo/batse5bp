batse5bp
=========

The batse5bp package provides persistant local access to BATSE GRB data hosted
at the CGRO SSC.  It accesses both the BATSE GRB catalogs, and a significant
amount of the detailed GRB data (e.g., time series and spectral data).

The last official BATSE GRB catalog was the 4B catalog.  The "current" catalog
hosted at the SSC includes many GRBs observed subsequent to the 4B catalog.  A
5B catalog, including data spanning to the end of the CGRO mission, was
planned but never completed (in particular, the sky exposure and efficiency
tables were never computed).  The "5bp" part of the package name denotes "5B,
preliminary."

The full "5Bp" GRB catalog is fetched, parsed, and stored locally when the
package is first used.  Detailed GRB data, which is much more voluminous, is
fetched transparently only when the user first attempts to access it; it is
cached locally for subsequent access.

Cached data is stored in a local file hierarchy with a structure resembling
that of the SSC archive.  Files and folders may be deleted within the local
filesystem when no longer needed; they will be re-fetched as necessary.  There
is no monolithic database file.

Documentation for the module is available online at:

http://inference.astro.cornell.edu/batse/

The documentation is build using the Sphinx documentation generator; the
Sphinx source is included in the batse5bp repository so the documentation
may be generated locally if needed.

NASA has supported development of the batse5bp package via its Research
Opportunities in Space and Earth Sciences (ROSES) programs.  Most of the
development has been supported by the ROSES Astrophysical Data Analysis Program
(ADAP) via grant NNX09AD03G.  Development of some of the time series data
capability has been additionally provided by the ROSES Applied Information
Systems Research Program (AISRP) via grant NNX09AK60G.
