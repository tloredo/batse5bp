#!/usr/bin/env python

# To install the package to a directory other than site-packages, use:
#
#     python setup.py install --install-platlib <path>

def configuration(parent_package='',top_path=''):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('batse5bp', parent_package, top_path)

    # Include subpackages:
    # config.add_subpackage('XXX')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        # name = 'batse5bp',  # use if name not given in configuration
        version = '0.1',
        packages = ['batse5bp'],
        install_requires=['numpy'],
        maintainer = "Tom Loredo",
        maintainer_email = "loredo@astro.cornell.edu",
        description = "Persistant access to BATSE GRB data hosted at the CGRO SSC",
        url = "http://www.notyetset.org/",
        license = 'LICENSE.txt',  # BSD license
        configuration=configuration)
