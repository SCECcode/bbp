#!/usr/bin/env python
"""
Copyright 2010-2020 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Correlation Matrix Generator for SDSU Broadband GroundMotion
Simulation Module on SCEC BBP (Python version)
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import numpy as np
import scipy as sp
from scipy.io import loadmat
from scipy import interpolate

# Import Broadband modules
from station_list import StationList

def get_number_of_points(input_file):
    """
    Reads input_file and returns the number of points in the file
    """
    npts = 0

    bbp_file = open(input_file, 'r')
    for line in bbp_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("%"):
            continue
        npts = npts + 1
    bbp_file.close()

    return npts

def get_dt(input_file):
    """
    Reads input_file and returns the dt used in the BBP file
    """
    val1 = None
    val2 = None

    ifile = open(input_file)
    for line in ifile:
        line = line.strip()
        # Skip blank lines
        if not line:
            continue
        # Skip comments
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        pieces = [float(piece) for piece in pieces]
        if val1 is None:
            val1 = pieces[0]
            continue
        if val2 is None:
            val2 = pieces[0]
            break
    ifile.close()

    # Quit if cannot figure out dt
    if val1 is None or val2 is None:
        print("[ERROR]: Cannot determine dt from file!")
        return None

    # Return dt
    return val2 - val1

def load_correlation_matrices(sdsu_data_dir):
    """
    Loads the .mat files with:
    - empirical correlation vectors and correlation matrix
    - empirical frequency vector and correlation matrices(3)
    """
    coeff_file1 = os.path.join(sdsu_data_dir, "empirical_coeff.mat")
    coeff_file2 = os.path.join(sdsu_data_dir, "B_EAS_2020.mat")
    data = {}

    data1 = loadmat(coeff_file1)
    data['emp_f'] = data1['emp_f'][0]
    data['emp_coeff'] = data1['emp_coeff']

    data2 = loadmat(coeff_file2)
    data['F_B'] = data2['F_B'][0]
    data['b'] = data2['b']
    data['b1'] = data2['b1']
    data['b2'] = data2['b2']
    data['b3'] = data2['b3']
    data['B'] = data2['B'][0]

    return data

def nearest_spd(input_matrix):
    """
    Calculates the nearest Symmetric Positive Definite Matrix
    """
    [r, c] = input_matrix.shape
    if r != c:
        print("ERROR: Input matrix should be a square matrix!")
        sys.exit(-1)

    # symmetrize A into B
    B = (input_matrix + input_matrix.T) / 2

    # Compute the symmetric polar factor of B. Call it H.
    # Clearly H is itself SPD
    U, s_diagonal, VH = np.linalg.svd(B)
    S = np.zeros((r, c))
    np.fill_diagonal(S, s_diagonal)
    V = VH.T
    H = V @ S @ V.T

    # get Ahat in the above formula
    Ahat = (B + H) / 2

    # Ensure symmetry
    Ahat = (Ahat + Ahat.T) / 2

    # Test that Ahat is in fact PD. if it is not so, then tweak it just a bit
    p = 1
    k = 0
    while p != 0:
        R, p = sp.linalg.lapack.dpotrf(Ahat)
        k = k + 1
        if p != 0:
            # Ahat failed the chol test, let's tweak it a bit
            my_eig, _ = np.linalg.eig(Ahat)
            min_eig = np.real(min(my_eig))
            Ahat = Ahat + (-min_eig * k**2 + np.spacing(min_eig)) * np.eye(r)

    return Ahat

def create_interfrequency_correlation(data):
    """
    For interfrequency correlation
    """
    print("Generating interfrequency correlation matrix...")

    #  Maximum and minimum empirical frequency of correlation implementation
    fcorrmax_inf = max(data['emp_f'])
    fcorrmin_inf = min(data['emp_f'])

    # Compute corresponding index of minimum empirical frequency and its range
    nks_inf = math.ceil(fcorrmin_inf / data['df'])
    # index in f that is less than empirical frequency value, should start
    # computation at nks+1 (f vector in BBcode starts at 0 Hz)
    nk_inf = math.floor(fcorrmax_inf / data['df']) + 1 - nks_inf
    # number of index in f that is within the empirical frequency value range

    # Generate interpolation mesh grid
    ff_inf = [a * data['df'] for a in range(nks_inf, nks_inf + nk_inf)]
    ffx_inf, ffy_inf = np.meshgrid(ff_inf, ff_inf)
    empfx_inf, empfy_inf = np.meshgrid(data['emp_f'], data['emp_f'])

    # Interpolation
    ip = interpolate.interp2d(empfx_inf[0], empfy_inf.T[0], data['emp_coeff'])
    Cint_inf =  ip(ffx_inf[0], ffy_inf.T[0])

    # Change diagonal elements to 1. Changing diagonal elements to 1 should
    # be able to make Cint positive definite in practical,
    # if not, use the following optional step.
    for i in range(0, len(ff_inf)):
        Cint_inf[i][i] = 1.0

    # Cholesky decomposition, keep the lower triangular matrix
    Kinf = np.linalg.cholesky(Cint_inf)

    return nk_inf, nks_inf, Kinf

def create_spatial_correlation(data):
    """
    For spatial correlation
    """
    print("Generating spatial correlation matrices...")

    # Maximum and minimum empirical frequency of correlation implementation
    fcorrmax_sp = max(data['F_B'])
    fcorrmin_sp = min(data['F_B'])

    # Compute corresponding index of minimum empirical frequency and its range

    # index in f that is less than empirical frequency value, should start
    # computation at nks+1 (f vector in BBcode starts at 0 Hz)
    nks_sp = math.ceil(fcorrmin_sp / data['df'])

    # Number of index in f that is within the empirical frequency value range
    nk_sp = math.floor(fcorrmax_sp / data['df']) + 1 - nks_sp

    # Generate interpolation mesh grid
    ff_sp = [a * data['df'] for a in range(nks_sp, nks_sp + nk_sp)]
    ffx_sp, ffy_sp = np.meshgrid(ff_sp, ff_sp)
    empfx_sp, empfy_sp = np.meshgrid(data['F_B'], data['F_B'])

    # Interpolation
    sp1 = interpolate.interp2d(empfx_sp[0], empfy_sp.T[0], data['b1'])
    sp2 = interpolate.interp2d(empfx_sp[0], empfy_sp.T[0], data['b2'])
    sp3 = interpolate.interp2d(empfx_sp[0], empfy_sp.T[0], data['b3'])
    Bint_sp1 = sp1(ffx_sp[0], ffy_sp.T[0])
    Bint_sp2 = sp2(ffx_sp[0], ffy_sp.T[0])
    Bint_sp3 = sp3(ffx_sp[0], ffy_sp.T[0])

    # Make sure Bint is positive definite
    Bint_sp1 = nearest_spd(Bint_sp1)
    Bint_sp2 = nearest_spd(Bint_sp2)
    Bint_sp3 = nearest_spd(Bint_sp3)

    # Cholesky decomposition, keep the lower triangular matrix
    Ksp1 = np.linalg.cholesky(Bint_sp1)
    Ksp2 = np.linalg.cholesky(Bint_sp2)
    Ksp3 = np.linalg.cholesky(Bint_sp3)

    return nk_sp, nks_sp, Ksp1, Ksp2, Ksp3

def generate_matrices(sdsu_data_dir, lf_seis_dir,
                      station_file, output_dir,
                      sim_id):
    """
    Generates the matrices
    """
    # Load correlation matrices
    data = load_correlation_matrices(sdsu_data_dir)

    # Figure out how long LF seismograms are
    slo = StationList(station_file)
    station_list = slo.getStationList()
    station_id = station_list[0].scode
    lf_file = os.path.join(lf_seis_dir, "%d.%s-lf.bbp" % (sim_id, station_id))
    ts_dt = get_dt(lf_file)
    wts_npts = get_number_of_points(lf_file)

    # Compute input LF time-series length, this is how BBcode
    # compute the LF timelength (starts at 0 second)
    wts_len = (wts_npts - 1) * ts_dt

    # Compute frequency step (df) used for correlation matrix elements

    # Fixed parameters
    tmp_npts = 32768
    tmp_lf_len = 102.3750

    # Compute total number of time point used in BB computation
    n = math.ceil(wts_len / tmp_lf_len)
    v_npts = (tmp_npts - 1) * n + 1
    npts = v_npts
    exponent = math.log(npts) / math.log(2.0)
    if exponent != np.fix(exponent):
        npts = 2**math.ceil(exponent)

    # Compute dt and df
    dt = tmp_lf_len / (tmp_npts - 1)
    data['df'] = 1 / (npts * dt)

    # For interfrequency correlation
    nk_inf, nks_inf, Kinf = create_interfrequency_correlation(data)

    # For spatial correlation
    nk_sp, nks_sp, Ksp1, Ksp2, Ksp3 = create_spatial_correlation(data)

    #if np.isnan(np.sum(Ksp1)):
    #    print("NaN in Ksp1!")
    #if np.isnan(np.sum(Ksp2)):
    #    print("NaN in Ksp2!")
    #if np.isnan(np.sum(Ksp3)):
    #    print("NaN in Ksp3!")
    #if np.isnan(np.sum(Kinf)):
    #    print("NaN in Kinf!")

    # Write output files
    kinf_fn = os.path.join(output_dir, "Kinf.bin")
    ksp1_fn = os.path.join(output_dir, "Ksp1.bin")
    ksp2_fn = os.path.join(output_dir, "Ksp2.bin")
    ksp3_fn = os.path.join(output_dir, "Ksp3.bin")

    kinf_fh = open(kinf_fn, 'wb')
    kinf_fh.write(nk_inf.to_bytes(4, 'little'))
    kinf_fh.write(nks_inf.to_bytes(4, 'little'))
    Kinf_32 = np.array(Kinf, dtype='float32')
    my_bytes = Kinf_32.tobytes('F')
    kinf_fh.write(my_bytes)
    kinf_fh.close()

    ksp1_fh = open(ksp1_fn, 'wb')
    ksp1_fh.write(nk_sp.to_bytes(4, 'little'))
    ksp1_fh.write(nks_sp.to_bytes(4, 'little'))
    Ksp1_32 = np.array(Ksp1, dtype='float32')
    my_bytes = Ksp1_32.tobytes('F')
    ksp1_fh.write(my_bytes)
    ksp1_fh.close()

    ksp2_fh = open(ksp2_fn, 'wb')
    Ksp2_32 = np.array(Ksp2, dtype='float32')
    my_bytes = Ksp2_32.tobytes('F')
    ksp2_fh.write(my_bytes)
    ksp2_fh.close()

    ksp3_fh = open(ksp3_fn, 'wb')
    Ksp3_32 = np.array(Ksp3, dtype='float32')
    my_bytes = Ksp3_32.tobytes('F')
    ksp3_fh.write(my_bytes)
    ksp3_fh.close()

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: %s sdsu_data_dir lf_seis_dir station_file output_dir sim_id" %
              (sys.argv[0]))
        sys.exit(-1)

    sdsu_data_dir = sys.argv[1]
    lf_seis_dir = sys.argv[2]
    station_file = sys.argv[3]
    output_dir = sys.argv[4]
    sim_id = int(sys.argv[5])

    generate_matrices(sdsu_data_dir, lf_seis_dir,
                      station_file, output_dir,
                      sim_id)
