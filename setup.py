import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='MAZSE',
    version='0.1.0',
    description='Multiscale Zeolite Atomic simulation Environment (MAZE)',
    url='https://github.com/kul-group/MAZE',
    author='Dexter Antonio',
    author_email='dexter.d.antonio@gmail.com',
    license="MIT",
    packages=['maze'],
    include_package_data=True,
    install_requires=['ase', 'numpy', 'typing'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: MIT',
        'Operating System :: MacOS, Windows, Linux',
        'Programming Language :: Python :: 3.8',
    ],

)
