#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Handles option file reading
"""
from __future__ import division, print_function

class OptFile(object):
    """
    This class is used to open and read an option file. It stores the
    information locally and returns each option individually, after each
    call to its get_next_option method.
    """

    def __init__(self, filename):
        """
        This function reads the option file pointed by filename.
        """
        self.filename = filename
        input_file = open(self.filename, 'r')
        self.data = input_file.readlines()
        input_file.close()
        self.line_no = 0

    def get_next_option(self):
        """
        This function returns the next option we read from the option file.
        """
        line = self.data[self.line_no]
        self.line_no += 1
        while line.startswith('#') or line.startswith('%'):
            line = self.data[self.line_no]
            self.line_no += 1
        comment_mark = line.find('#')
        if comment_mark < 0:
            comment_mark = line.find('%')
        if comment_mark >= 0:
            return line[0:comment_mark].strip()
        else:
            return line.strip()
