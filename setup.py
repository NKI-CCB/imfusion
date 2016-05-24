import sys
import setuptools

import versioneer


# General requirements.
INSTALL_REQUIRES = ['future', 'pandas', 'numexpr', 'pysam',
                    'toolz', 'scipy', 'seaborn', 'pyfaidx',
                    'scipy', 'intervaltree']

EXTRAS_REQUIRE = {
    'de_single': ['rpy2'],
    'dev': ['sphinx', 'pytest', 'pytest-mock',
            'pytest-datafiles', 'pytest-cov',
            'pytest-helpers-namespace']
}


# Check setuptools version, as recommended by:
# https://hynek.me/articles/conditional-python-dependencies/.
if int(setuptools.__version__.split('.', 1)[0]) < 18:
    assert 'bdist_wheel' not in sys.argv

    # Add pathlib for Pythons before 3.4.
    if sys.version_info[0:2] < (3, 4):
        INSTALL_REQUIRES.append('pathlib2')
else:
    EXTRAS_REQUIRE[":python_version<'3.4'"] = ['pathlib2']


setuptools.setup(
    name='im-fusion',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url='https://github.com/jrderuiter/im-fusion',
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    description=('Tool for identifying transposon insertions in '
                 'Insertional Mutagenesis screens from gene-transposon '
                 'fusions using single- and paired-end RNA-sequencing data.'),
    license='MIT',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    entry_points={'console_scripts': [
        'im-fusion = imfusion.main:main',
    ]},
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    zip_safe=True,
    classifiers=[]
)
