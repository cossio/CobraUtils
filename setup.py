'''
pip install
'''

import glob
from setuptools import setup


setup(name='CobraUtils', version='0.1.0', license='MIT',
      url='https://github.com/cossio/CobraUtils',
      author='Jorge Fernandez-de-Cossio-Diaz',
      author_email='j.cossio.diaz@gmail.com',
      description='Utilities for the cobrapy package.',
      packages=['CobraUtils'],
      scripts=glob.glob("bin/*.py") + glob.glob("bin/*.sh"),
      long_description=open('README.md').read(),
      zip_safe=False,
      setup_requires=['cobra'],
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True)
