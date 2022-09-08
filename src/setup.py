
from setuptools import setup, find_packages, Extension

setup(
    name='spkmodule',
    version='0.1.0',
    author='TalAndMaya',
    author_email='tal@none.com',
    description='Final project wrapper',
    install_requires=['invoke'],
    packages=find_packages(where='.', exclude=()),
    license='GPL-2',
    classifiers=[
        'development Status :: prod',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'spkmodule',
            ['spkmeansmodule.c'],
        )
    ]
)
