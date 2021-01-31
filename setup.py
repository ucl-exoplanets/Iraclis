
import os
import glob


from setuptools import setup

package = 'iraclis'
version = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__version__.txt')).read()
author = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__author__.txt')).read()
author_email = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__author_email__.txt')).read()
description = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__description__.txt')).read()
url = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__url__.txt')).read()

install_requires = ['docopt', 'pylightcurve>=4.0.0', 'numpy>=1.19.2', 'matplotlib>=3.3.2', 'scipy>=1.5.2',
                    'emcee>=3.0.2', 'sklearn', 'astropy>=4.2']
entry_point = '__main__:console'

os.chdir(os.path.abspath(os.path.dirname(__file__)))

subdirs_to_include = []
for x in os.walk(package):
    if os.path.isdir(x[0]):
        if x[0] != package:
            subdirs_to_include.append(x[0])

files_to_include = []
for x in glob.glob(os.path.join(package, '*')):
    if os.path.isfile(x):
        if x.split('.')[-1] not in ['py']:
            files_to_include.append(os.path.join(package, os.path.split(x)[1]))

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

setup(
    name=package,
    version=version,
    description=description,
    long_description='Visit {0} for more details.'.format(url),
    url=url,
    author=author,
    author_email=author_email,
    license='Creative Commons Attribution 4.0 International License',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X',
                 'Programming Language :: Python :: 3.7',
                 ],
    entry_points={'console_scripts': ['{0} = {0}.{1}'.format(package, entry_point)]},
    packages=[package],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
    setup_requires=["pytest-runner"],
    tests_require=['pytest'],
)
