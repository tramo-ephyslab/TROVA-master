from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
import setuptools
import subprocess
import sysconfig

with open("README.md", "r") as fh:
    long_description = fh.read()

ext_suffix = sysconfig.get_config_var('EXT_SUFFIX') or '.so'
trova_so = os.path.join('functions' + ext_suffix)


try:
	print(subprocess.check_output(
		"cd trova/; "
		"sh buid_trova_so.sh", shell=True))
except subprocess.CalledProcessError as e:
	print(e.output)
	print("Problems compiling the function module.  "
			"Will continue using a slower fallback...")
else:

	print()

	print("-----------------------------------------------------------------------")
	print("FORTRAN functions extension has been created as {}".format(trova_so))
	print("-----------------------------------------------------------------------")


setuptools.setup(
    name="trova",
    version="1.1",
    author="José C. Fernández-Alvarez & Albenis Pérez-Alarcón",
    author_email="jose.carlos.fernandez.alvarez@uvigo.es & albenis.pérez.alarcon@uvigo.es",
    description="TROVA is a Python software for studying moisture sources and sinks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tramo-ephyslab/TROVA-master",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy","mpi4py","time","struct","datetime","netCDF4","scipy","functools","pathlib","gzip","shutil","imp", "matplotlib", "cartopy", "math"],
    include_package_data=True,
    package_data={"":['*.so','VERSION', '*.f90',"_version.py","LAST_UPDATE"]},

)
