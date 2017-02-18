#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

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
