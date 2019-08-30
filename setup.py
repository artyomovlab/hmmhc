from setuptools import setup

setup(
    name='hmmhc',
    version='0.1',
    description='Hidden Markov model-based MHC II binding predictor',
    author='Maxim Artyomov and Ilya Kizhvatov',
    author_email='ilya.kizhvatov@wustl.edu',
    url='https://github.com/artyomovlab/hmmhc',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
    packages=['hmmhc'],
    install_requires=[
        'numpy',
        'pandas'
    ],
    entry_points = {
        'console_scripts': ['hmmhc-predict=hmmhc.cmdline:main'],
    },
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False
)
