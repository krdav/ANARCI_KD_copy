#!/usr/bin/env python

from distutils.core import setup
setup(name='anarci',
      version='1.1',
      description='Antibody Numbering and Receptor ClassIfication',
      author='James Dunbar',
      author_email='dunbar@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/ANARCI',
      packages=['anarci', 'anarci.Bio', 'anarci.Bio.Align', 
      'anarci.Bio.Alphabet', 'anarci.Bio.Data', 'anarci.Bio._py3k',
      'anarci.Bio.SearchIO','anarci.Bio.SearchIO._model', 'anarci.Bio.SearchIO.HmmerIO'  ],
      package_dir={'anarci': 'lib/python/anarci', 
                   'anarci.Bio.SearchIO': 'lib/python/anarci/Bio/SearchIO', 
                   'anarci.Bio.SearchIO.HmmerIO': 'lib/python/anarci/Bio/SearchIO/HmmerIO', 
                   'anarci.Bio.Alphabet': 'lib/python/anarci/Bio/Alphabet', 
                   'anarci.Bio._py3k': 'lib/python/anarci/Bio/_py3k', 
                   'anarci.Bio.Data': 'lib/python/anarci/Bio/Data', 
                   'anarci.Bio': 'lib/python/anarci/Bio', 
                   'anarci.Bio.SearchIO._model': 'lib/python/anarci/Bio/SearchIO/_model', 
                   'anarci.Bio.Align': 'lib/python/anarci/Bio/Align'},
      package_data={'anarci': ['dat/HMMs/ALL.hmm',
                               'dat/HMMs/ALL.hmm.h3f',
                               'dat/HMMs/ALL.hmm.h3i',
                               'dat/HMMs/ALL.hmm.h3m',
                               'dat/HMMs/ALL.hmm.h3p']},
      scripts=['bin/ANARCI'],
      license="GPLv3"
     )
