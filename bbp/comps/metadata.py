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

This program defines a metadata class for managing metadata
files used in workflows. This class will be used by workflow
components.
"""
from __future__ import division, print_function

class Metadata(object):

    def __init__(self, base_prog_name):
        """
        Programs shoudl pass in a string that identifies. This string is
        then pre-pended to all attributes.
        """
        self.mdict = {}
        self.base = base_prog_name
        self.mdict["workflow_name"] = self.base
        self.seqnum = 0

    def add(self, attrib, value):
        """
        This method will input a metadatda attribute and a metadata value.
        These will get prepended to the metadata dictioinary.

        """
        self.seqnum = self.seqnum + 1
        attr = "%s.%05d_%s" % (self.base, self.seqnum, attrib)
        self.mdict[attr] = value

    def create_file(self, a_metadata):
        f = open(a_metadata, "w")
        f.write("#Starting a Metadata File for the Broadband Platform\n")
        for i, v in self.mdict.items():
            line = "%s=%s\n" % (i, v)
            f.write(line)
        f.close()

    def show_metadata(self):
        for i, v in self.mdict.items():
            print("%s : %s\n" % (i, v))

    def write_to_file(self, a_metadata):
        """
        This method will write the xml to a metadata file as specified
        on the command line
        """
        f = open(a_metadata, "a")
        for i, x in self.mdict.items():
            attrib = i
            value = x
            val = "%s=%s\n" % (attrib, value)
            f.write(val)
        f.close

if __name__ == "__main__":
    meta = Metadata("jbsim3d")
    meta.add("author", "Philip Maechling")
    meta.show_metadata()
