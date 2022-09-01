from setuptools import setup, find_packages
import glob
import os
import pkg_resources
import pathlib

from oat import __version__, _program

# glob_fix command thanks to https://stackoverflow.com/a/64789489/15793575
def glob_fix(package_name, glob):
    # this assumes setup.py lives in the folder that contains the package
    package_path = pathlib.Path(f'./{package_name}').resolve()
    return [str(path.relative_to(package_path))
            for path in package_path.glob(glob)]

setup(name='oat',
      version=__version__,
      packages=['oat', 'oat.scripts'],
#      packages=find_packages(),
      scripts=['oat/scripts/analysis_module.smk','oat/scripts/alternate_analysis.smk','oat/scripts/get_latest_tag.sh'
                ],
      package_data={"oat":[*glob_fix("oat", "protocols/**/*"),
                            *glob_fix("oat", "references/**/*"),
                            *glob_fix("oat", "envs/**/*")]},
     # package_data={"oat":["protocols/**/*","references/**/*"]},
      install_requires=[
            'pandas>=1.0.1',
        ],
      description='ont-analysis-toolkit from SAVID: Snakemake edition',
      url='https://github.com/charlesfoster/ont-analysis-toolkit',
      author='Dr Charles Foster',
      author_email='charles.foster@unsw.edu.au',
      entry_points="""
      [console_scripts]
      {program} = oat.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
