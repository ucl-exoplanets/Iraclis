from setuptools import setup
import codecs
import os
import glob

name = 'iraclis'
description = 'Analysis pipeline for HST/WFC3 spectroscopic observations of exoplanet transits and eclipses'
url = 'https://github.com/ucl-exoplanets/Iraclis'
install_requires = ['docopt', 'pylightcurve>=3.0.5', 'numpy>=1.19.2', 'matplotlib>=3.3.2', 'scipy>=1.5.2',
                    'emcee>=3.0.2', 'sklearn', 'astropy>=4.2']
entry_point = '__main__:console'

os.chdir(os.path.abspath(os.path.dirname(__file__)))

subdirs_to_include = []
for x in os.walk(name):
    if os.path.isdir(x[0]):
        if x[0] != name:
            subdirs_to_include.append(x[0])

files_to_include = []
for x in glob.glob(os.path.join(name, '*')):
    if os.path.isfile(x):
        if x.split('.')[-1] not in ['py']:
            files_to_include.append(os.path.join(name, os.path.split(x)[1]))

files_to_include.append('README.md')
files_to_include.append('LICENSE')
files_to_include.append('readme.md')
files_to_include.append('licence')

w = open('MANIFEST.in', 'w')
for i in subdirs_to_include:
    w.write('include ' + os.path.join(i, '*') + ' \n')

for i in files_to_include:
    w.write('include ' + i + ' \n')

w.close()

try:
    with codecs.open('README.md', encoding='utf-8') as f:
        long_description = f.read()
except:
    with codecs.open('readme.md', encoding='utf-8') as f:
        long_description = f.read()

version = ' '
for i in open(os.path.join(name, '__init__.py')):
    if len(i.split('__version__')) > 1:
        version = i.split()[-1][1:-1]

setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    url=url,
    author='Angelos Tsiaras',
    author_email='aggelostsiaras@gmail.com',
    license='Creative Commons Attribution 4.0 International License',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X',
                 'Programming Language :: Python :: 3.7',
                 ],
    entry_points={'console_scripts': ['{0} = {0}.{1}'.format(name, entry_point)]},
    packages=[name],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)
