from setuptools import setup

setup(
    name='MAZE',
    version='0.1.0',
    description='Multiscale Zeolite Simulation Enviorment (MAZE)',
    url='https://github.com/shuds13/pyexample',
    author='Dexter Antonio',
    author_email='dexter.d.antonio@gmail.com',
    license='Private',
    packages=['pyexample'],
    install_requires=['ase',
                      'numpy'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'Operating System :: MacOS or Windows',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.9',
    ],
)