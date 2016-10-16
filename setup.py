#!/usr/bin/env python
# coding: utf-8
from setuptools import setup, find_packages
from isoclique import __author__, __version__, __license__
 
setup(
        name             = 'isoclique',
        version          = __version__,
        description      = 'Implementation of Isolated Clique Enumeration Algorithm proposed in "Linear-Time Enumeration of Isolated Cliques". See http://link.springer.com/chapter/10.1007%2F11561071_13 for algorithm details.',
        license          = __license__,
        author           = __author__,
        author_email     = 'ahasimoto@mm.media.kyoto-u.ac.jp',
        url              = 'https://github.com/AtsushiHashimoto/isolated_clique_enumeration.git',
        keywords         = 'isolated cliques',
        packages         = find_packages(),
        install_requires = ['numpy'],
        )
