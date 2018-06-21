# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

import os
import shutil
import sys
try:
    from urllib.parse import urljoin
except ImportError:
    from urlparse import urljoin

import requests


def main(token, filepath):
    assert 'CPLEX_URL' in os.environ
    url = urljoin(os.environ['CPLEX_URL'], 'contents')
    sha_url = urljoin(os.environ['CPLEX_URL'], 'git/blobs/{}')
    headers = {
        'Authorization': 'token {}'.format(token),
        'Accept': 'application/vnd.github.v3.raw'
    }
    files_meta = requests.get(url, headers=headers).json()
    sha = [i for i in files_meta if i['path'] == filepath][0]['sha']
    response = requests.get(sha_url.format(sha), headers=headers, stream=True)
    with open(filepath, 'wb') as out_file:
        response.raw.decode_content = True
        shutil.copyfileobj(response.raw, out_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:\n{} <GitHub Token> <Archive Name>")
        sys.exit(2)
    main(*sys.argv[1:])
