"""
Provide access to DRM data stored in FITS files at COSSC.

Documentation of the compressed DRM format from the DRM FITS header
(for 1-based array indexing):

The detector response matrix is stored in a compressed format:        

Elements are listed in order of row by column (n_e_bin by n_e_chan), with
all elements equal to 0 at the top of the column ommited.

The row index of the first non-zero element for each column is given by
n_zeros.

Thus, for the ith column: insert n_zeros(i) - 1 0s, followed by the next
n_e_bin - n_zeros(i) + 1 elements of the drm (either drm_dir, drm_sct or
drm_sum).

The last column should exhaust the elements of the compressed drm.

Created 2012-10-16 by Tom Loredo
"""
from os.path import join, exists
import cPickle

from numpy import zeros, minimum

import pyfits


# This list of the binary table header fields was obtained from a DISCSC DRM
# file by examining the binary table header __dict__.
hdr_fields = [
    ('TTYPE1', 'DET_NUM', 'the active detector for one row'),
    ('TFORM1', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE2', 'MAT_TYPE', '0=Direct,1=Scattered, 2= Both, 3=Summed'),
    ('TFORM2', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE3', 'TIME', 'time drm is calculated for'),
    ('TFORM3', '1D', 'data format of the field: DOUBLE PRECISION'),
    ('TUNIT3', 'TJD', 'physical unit field'),
    ('TTYPE4', 'NUMEBINS', 'number of rows in PHT_EDGE (N_E_BINS+1)'),
    ('TFORM4', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE5', 'NUMCHAN', 'number of columns in E_EDGES (N_E_CHAN+1)'),
    ('TFORM5', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE6', 'NUMZERO', 'number of of zeroes top of each column'),
    ('TFORM6', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE7', 'DIRDRM', 'number of rows in the DIR drm'),
    ('TFORM7', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE8', 'SCTDRM', 'number of rows in the  SCT drm'),
    ('TFORM8', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE9', 'SUMDRM', 'number of rows in the  SUM drm'),
    ('TFORM9', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE10', 'PHT_EDGE', 'input bin energy edges'),
    ('TFORM10', 'PE(276)', 'data format of the field'),
    ('TUNIT10', 'keV', 'physical unit field'),
    ('TTYPE11', 'E_EDGES', 'channel energy edges'),
    ('TFORM11', 'PE(253)', 'data format of the field'),
    ('TUNIT11', 'keV', 'physical unit field'),
    ('TTYPE12', 'N_ZEROS', 'list of how many 0 elements in each column'),
    ('TFORM12', 'PI(252)', 'data format of the field'),
    ('TTYPE13', 'DRM_DIR', 'detector response matrix'),
    ('TFORM13', 'PE(69552)', 'data format of the field'),
    ('TUNIT13', 'cm**2', 'physical unit field'),
    ('TTYPE14', 'DRM_SCT', 'detector response matrix'),
    ('TFORM14', 'PE(69552)', 'data format of the field'),
    ('TUNIT14', 'cm**2', 'physical unit field'),
    ('TTYPE15', 'DRM_SUM', 'detector response matrix'),
    ('TFORM15', 'PE(69552)', 'data format of the field'),
    ('TUNIT15', 'cm**2', 'physical unit field'),
    ('TTYPE16', 'CAL_NAME', 'Name of energy calibration method'),
    ('TFORM16', '16A', 'data format of the field') ]


# Collect information for the TTYPE fields in a dict.
ttype_fields = {}
for attr, val, comment in hdr_fields:
    if attr[:5] == 'TTYPE':
        col = int(attr[5:]) - 1
        ttype_fields[col] = (val, comment)


class DRM_DISCSC:
    """
    Provide access to the detector response matrix (DRM) for a single BATSE
    detector's discriminator science (DISCSC) data.

    The DRM units according to the FITS header is area in cm^2.
    """

    fields = ttype_fields

    def __init__(self, row_data, n_ch, n_E):
        """
        Retrieve DRM info from a single row of parsed binary FITS table data.
        """ 
        for col, info in self.fields.items():
            setattr(self, info[0], row_data[col])

        self.det_num = self.DET_NUM

        # Discriminator bin boundaries (nominal energy loss units):
        self.ch_bins = self.E_EDGES  # n_ch+1 values

        # Incident photon energy "edges;" there should be n_E of these but
        # they are treated as if binned somehow....
        self.E_bins = self.PHT_EDGE

        # Unpack the zero-shortened format.
        self.start = self.N_ZEROS - 1  # N_ZEROS is 1-based index of 1st non-0 element
        self.drm = zeros((n_ch, n_E))
        run = 0
        for i in range(n_ch):
            num_nz = n_E - self.start[i]
            # Use the summed (direct + atmospheric scattering) matrix; typically
            # the direct and scattering entries are empty.
            self.drm[i,self.start[i]:] = self.DRM_SUM[run:run+num_nz]
            run += num_nz

        # Make views providing detector channel response vectors.
        self.crv = []
        for i in range(n_ch):
            self.crv.append(self.drm[i,:])


class DRMs_DISCSC:
    """
    Provide access to detector response matrix (DRM) data from a BATSE DISCSC DRM
    FITS file for all triggered detectors associated with a BATSE trigger; also
    calculate the summed response.
    """

    # Metadata attributes stored in the "-meta" file.
    meta_attrs = ['burst_id', 'file_id', 'n_ch', 'n_E', 'alpha', 'n_det']

    def __init__(self, grb):
        """
        Load DRM data for the DISCSC data associated with GRB instance `grb`.

        If the DRM data is not locally cached, it is fetched from the SSC,
        loaded, and cached for future use.  Otherwise, the cached copy is
        loaded.
        """
        self.grb = grb

        # Get DRMs for the triggered detectors.
        drm_path = 'discsc_drm_%i.pkl' % grb.trigger  # stores detector DRMs
        drm_path = join(self.grb.grb_dir, drm_path)
        meta_path = 'discsc_drm_%i-meta.pkl' % grb.trigger  # stores metadata
        meta_path = join(self.grb.grb_dir, meta_path)
        if exists(drm_path) and exists(meta_path):  # already parsed FITS
            ifile = open(meta_path, 'rb')
            metadata = cPickle.load(ifile)
            ifile.close()
            for attr in self.meta_attrs:
                setattr(self, attr, metadata[attr])
            ifile = open(drm_path, 'rb')
            self.detectors = cPickle.load(ifile)
            ifile.close()
            try:
                assert self.n_det == len(self.detectors)
            except:
                raise RuntimeError('Data mismatch in archive DRM files!')
        else:  # must parse FITS
            fits_name = 'discsc_drm_%i.fits.gz' % grb.trigger
            fits_path = grb.raw_cached_path(fits_name)
            self.load_from_fits(fits_path, drm_path, meta_path)

        # The full-instrument DRM sums the triggered detector DRMs.
        self.drm = zeros((self.n_ch, self.n_E))
        for i, d in enumerate(self.detectors):
            self.drm += d.drm
            # Keep track of lowest non-zero entry for each summed channel.
            if i == 0:
                self.start = d.start[:]
            else:
                self.start = minimum(self.start, d.start)

        # Make views providing detector channel response vectors.
        self.crv = []
        for i in range(self.n_ch):
            self.crv.append(self.drm[i,:])

        # Get energy node info for sum DRM from one of the detectors.
        # Discriminator bin boundaries (nominal energy loss units), n_ch+1 vals:
        self.ch_bins = self.detectors[0].ch_bins
        # Incident photon energy "edges;" there should be n_E of these but
        # they are treated as if binned somehow....
        self.E_bins = self.detectors[0].E_bins

    def load_from_fits(self, fits_path, drm_path, meta_path):
        """
        Read DRM data from a FITS file.  The path `fname` may be to a gzipped
        version.  Archive the data at the DRM file paths.
        """
        hdus = pyfits.open(fits_path)
        primary = hdus[0].header

        # Verify some key primary header elements.
        try:
            assert primary['filetype'] == 'BATSE_DRM'
            assert primary['det_mode'] == 'LAD'
        except:
            raise ValueError('Primary header is not as expected for DRM FITS file!')

        self.burst_id = primary['object']  # GRByymmdd designator
        self.file_id = primary['file-id']  # contains trigger #

        self.n_ch = primary['n_e_chan']  # number of discriminator channels
        # number of of incident photon energies (*not* bins):
        self.n_E = primary['n_e_bins']
        self.alpha = primary['alpha']  #  weighting across input bins for direct matrix

        # Read the DRM binary table extension, close the file, and gather
        # some basic info (field names and sizes).
        exten = hdus[1].header
        try:
            assert exten['exttype'] == 'BATSEDRM'
        except:
            raise ValueError('Table header is not as expected for DRM FITS file!')

        # These are handy for interactive exploring/checking of FITS data:
        # self.primary = primary
        # self.exten = exten
        # self.data = hdus[1].data

        # Each row contains the DRM for one of the triggered detectors.
        self. n_det = len(hdus[1].data)  # number of triggered detectors
        self.detectors = []
        for row in hdus[1].data:
            self.detectors.append(DRM_DISCSC(row, self.n_ch, self.n_E))
        hdus.close()

        # Store the data in a format for easier reloading.
        metadata = {}
        for attr in self.meta_attrs:
            metadata[attr] = getattr(self, attr)
        ofile = open(meta_path, 'wb')
        cPickle.dump(metadata, ofile)
        ofile.close()
        ofile = open(drm_path, 'wb')
        cPickle.dump(self.detectors, ofile, -1)  # use highest protocol
        ofile.close()

