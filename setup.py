from setuptools import setup, find_packages, Extension
import os
import shutil
import sys

# Get the path to the bin directory in the Conda environment
bin_dir = os.path.join(sys.prefix, 'bin')


setup(
    name="localhgt",
    version="1.0.0",
    author="wangshuai",
    author_email="wshuai294@gmail.com",
    description="Detect HGT from microbiome.",
    url="https://github.com/deepomicslab/LocalHGT", 
    packages=find_packages(),

    # ext_modules = [
    #     Extension(
    #         name = 'extract_ref',
    #         sources = [
    #             'src/extract_ref_normal_peak.cpp',
    #         ],
    #         include_dirs = ['src'],
    #     )
    # ],
    # data_files = [ ('bin', ['bin/extract_ref']) ],

    # data_files=[('/home/wangshuai/miniconda3/envs/hgt/bin', ['scripts/extract_ref'])],
    data_files=[(bin_dir, ['scripts/extract_ref'])],

    entry_points={
    'console_scripts': [
        'localhgt = localhgt:main',
    ],
    },

    scripts=['scripts/infer_HGT_breakpoint.py','scripts/infer_HGT_event.py', 'scripts/localhgt.py', 'scripts/accurate_bkp.py', 'scripts/extractSplitReads_BwaMem.py', 'scripts/get_bed_file.py', 'scripts/get_raw_bkp.py', 'scripts/pipeline.sh', 'scripts/remove_repeat.py']
)
