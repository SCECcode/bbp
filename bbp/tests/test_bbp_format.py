#! /usr/bin/env python
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

Purpose: These are acceptance tests for the PEER seismogram converter
Author: Philip Maechling
Revision Date: 12.7.27
"""
from __future__ import division, print_function

# Import Python modules
import os
import unittest

# Import Broadband modules
import seqnum
import bband_utils
import bbp_formatter
from install_cfg import InstallCfg

class TestBBPFormat(unittest.TestCase):
    """
    Unit tests for the bbp_formatter module
    """

    def setUp(self):
        """
        Configures the environment for the tests
        """
        self.install = InstallCfg.getInstance()
        self.sim_id = int(seqnum.get_seq_num())

        # Make sure all directories exist
        self.indir = os.path.join(self.install.A_IN_DATA_DIR,
                                  str(self.sim_id))
        self.tmpdir = os.path.join(self.install.A_TMP_DATA_DIR,
                                   str(self.sim_id))
        self.outdir = os.path.join(self.install.A_OUT_DATA_DIR,
                                   str(self.sim_id))
        self.logdir = os.path.join(self.install.A_OUT_LOG_DIR,
                                   str(self.sim_id))
        bband_utils.mkdirs([self.indir, self.tmpdir, self.outdir, self.logdir],
                           print_cmd=False)

    def test_bbp2peer(self):
        """
        Test for the bbp2peer converter
        """
        ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucb")
        in_bbp_file = os.path.join(ref_dir, "station.acc.bbp")
        peer_n_out = os.path.join(self.outdir, "output_peer_n.acc")
        peer_e_out = os.path.join(self.outdir, "output_peer_e.acc")
        peer_z_out = os.path.join(self.outdir, "output_peer_z.acc")
        for out_file in [peer_n_out, peer_e_out, peer_z_out]:
            try:
                os.remove(out_file)
            except:
                pass

        #
        # Simplest reference filename.
        #
        peer_n_out_ref = os.path.join(ref_dir, "station.peer_n.acc")
        peer_e_out_ref = os.path.join(ref_dir, "station.peer_e.acc")
        peer_z_out_ref = os.path.join(ref_dir, "station.peer_z.acc")

        bbp_formatter.bbp2peer(in_bbp_file, peer_n_out,
                               peer_e_out, peer_z_out)
        print("created PEER-format files: %s %s %s" %
              (peer_n_out, peer_e_out, peer_z_out))

        res_file = open(peer_n_out, 'r')
        lines = res_file.readlines()
        res_file.close()

        ref_file = open(peer_n_out_ref, 'r')
        rlines = ref_file.readlines()
        ref_file.close()
        #
        # Pass test if bbp file has more than 400 lines.
        # Only confirms that the bbp file is non-zero length
        #
        self.assertTrue(len(lines) == len(rlines))
        # convert to bbp format
        # input to respect
        # retrieve amplitudes from respect
        # compare to amplitudes calculated by PEER
        return 0

    def test_peer2bbp(self):
        """
        Test for the peer2bbp converter
        """
        ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "ucb")
        ## Retrieve PEER seismogram files
        in_n_file = os.path.join(ref_dir, "station.peer_n.acc")
        in_e_file = os.path.join(ref_dir, "station.peer_e.acc")
        in_z_file = os.path.join(ref_dir, "station.peer_z.acc")

        # Output file
        peer_out_bbp = os.path.join(self.outdir, "station.output.bbp")
        try:
            os.remove(peer_out_bbp)
        except:
            pass
        # Reference bbp filename
        peer_out_bbp_ref = os.path.join(ref_dir, "station.acc.bbp")

        bbp_formatter.peer2bbp(in_n_file, in_e_file,
                               in_z_file, peer_out_bbp)
        print("created bbp file: %s" % (peer_out_bbp))

        res_file = open(peer_out_bbp, 'r')
        lines = res_file.readlines()
        res_file.close()

        ref_file = open(peer_out_bbp_ref, 'r')
        rlines = ref_file.readlines()
        ref_file.close()
        #
        # Pass test if bbp file has more than 400 lines.
        # Only confirms that the bbp file is non-zero length
        #
        self.assertTrue(len(lines) == len(rlines))
        # convert to bbp format
        # input to respect
        # retrieve amplitudes from respect
        # compare to amplitudes calculated by PEER
        return 0

    def test_exsim2bbp(self):
        """
        Test for the exsim2bbp converter
        """
        ref_dir = os.path.join(self.install.A_TEST_REF_DIR, "uwo")
        ## Retrieve ExSIM seismogram files
        in_n_file = os.path.join(ref_dir, "LOMAPS001iter001.acc")
        in_e_file = os.path.join(ref_dir, "LOMAPS001iter002.acc")
        in_z_file = os.path.join(ref_dir, "LOMAPS001iter003.acc")

        # Output file
        peer_out_bbp = os.path.join(self.outdir, "lomap-output.bbp")
        try:
            os.remove(peer_out_bbp)
        except:
            pass
        #
        # Reference bbp filename
        peer_out_bbp_ref = os.path.join(ref_dir, "LOMAPS001.acc.bbp")

        bbp_formatter.exsim2bbp(in_n_file, in_e_file, in_z_file, peer_out_bbp)
        print("created bbp file: %s" % (peer_out_bbp))

        res_file = open(peer_out_bbp, 'r')
        lines = res_file.readlines()
        res_file.close()

        ref_file = open(peer_out_bbp_ref, 'r')
        rlines = ref_file.readlines()
        ref_file.close()
        #
        # Pass test if bbp file has more than 400 lines.
        # Only confirms that the bbp file is non-zero length
        #
        self.assertTrue(len(lines) == len(rlines))
        # convert to bbp format
        # input to respect
        # retrieve amplitudes from respect
        # compare to amplitudes calculated by PEER
        return 0

if __name__ == '__main__':
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestBBPFormat)
    unittest.TextTestRunner(verbosity=2).run(SUITE)
