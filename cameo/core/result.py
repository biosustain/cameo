# -*- coding: utf-8 -*-
# Copyright 2015 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import getpass
import time
from datetime import datetime

from cameo import system_info


class MetaInformation(object):
    def __init__(self, *args, **kwargs):
        super(MetaInformation, self).__init__(*args, **kwargs)
        self._system_info = system_info
        self._responsible = getpass.getuser()
        self._timestamp = time.time()

    @property
    def system_info(self):
        return self._system_info

    @property
    def responsible(self):
        return self._responsible

    @property
    def timestamp(self):
        """doc string"""
        return self._timestamp

    @property
    def human_readable_timestamp(self):
        dt = datetime.fromtimestamp(self.timestamp)
        return dt.strftime('%Y-%m-%d %H:%M:%S.%f')


class Result(object):
    def __init__(self, *args, **kwargs):
        super(Result, self).__init__(*args, **kwargs)
        self._meta_information = MetaInformation()

    @property
    def meta_information(self):
        return self._meta_information

    @property
    def data_frame(self):
        raise NotImplementedError

    def plot(self, grid=None, width=None, height=None, title=None, *args, **kwargs):
        raise NotImplementedError
