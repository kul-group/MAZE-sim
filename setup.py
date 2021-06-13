import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='MAZE-sim',
    version='0.1.1',
    description='Multiscale Zeolite Atomic simulation Environment (MAZE)',
    long_description="This project aims to extend the Atomic Simulation Environment (ASE) "
                     "to more naturally represent the properties of zeolites and facilitate "
                     "the calculations required to determine their properties. ",

    long_description_content_type='text/x-rst',
    url='https://github.com/kul-group/MAZE',
    author='Dexter Antonio',
    author_email='dexter.d.antonio@gmail.com',
    license="MIT",
    packages=['maze'],
    include_package_data=True,
    install_requires=['ase', 'numpy', 'typing', 'packaging', 'scikit-learn'],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='ase zeolites simulations',
    project_urls={
        'Documentation': 'https://kul-group.github.io/MAZE-sim/build/html/index.html',
        'Source': 'https://github.com/kul-group/MAZE-sim/',
        'Tracker': 'https://github.com/kul-group/MAZE-sim/issues',
    },
    python_requires='>=3'
)
