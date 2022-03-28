#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
