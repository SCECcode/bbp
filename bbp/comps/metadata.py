#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
