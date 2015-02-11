__all__ = ['bigg2mnx', 'mnx2bigg']

import os
import pickle
import gzip

import cameo

with open(os.path.join(cameo._cameo_data_path, 'metanetx.pickle')) as f:
    _METANETX = pickle.load(f)

bigg2mnx = _METANETX['bigg2mnx']
mnx2bigg = _METANETX['mnx2bigg']

with gzip.open(os.path.join(cameo._cameo_data_path, 'metanetx_chem_prop.pklz')) as f:
    chem_prop = pickle.load(f)

