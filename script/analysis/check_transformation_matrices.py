# ======================================================================
# copyright 2020. Triad National Security, LLC. All rights
# reserved. This program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
# is operated by Triad National Security, LLC for the U.S. Department
# of Energy/National Nuclear Security Administration. All rights in
# the program are reserved by Triad National Security, LLC, and the
# U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others
# acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# ======================================================================

# Authors: Oleg Korobkin (korobkin@lanl.gov)
# Purpose:
# Provides a check of whether a coordinate transformation of the metric
# from code coordinates to Kerr-Schild coordinates produces correct 
# metric, consistent with the closed form (as in e.g. Eq.(3)
# McKinney & Gammie 2004, https://arxiv.org/abs/astro-ph/0404512)
#
# Functions:
# - print_matrix
# - check_transformation_matrices
# 

from math import *
import numpy as np

def print_matrix(matrix,fmt="%19.11e",tostdout=True) -> str:
  """Pretty-prints a matrix to a string (optinally, to stdout)

  Parameters
  ----------
  matrix : numpy.array([N,M])
      matrix to print
  fmt : str
      C-style format of each element (default: "%19.11e")
  tostdout : bool
      output to stdout (default: true)

  Returns
  -------
  str
      formatted output string
  """
  N = matrix.shape[0]
  M = matrix.shape[1]
  s = "["
  for i in range(N):
    s+= "["
    for j in range(M):
      s+= (fmt % matrix[i,j])
      if j < M - 1: s += ", "
    s+= "]"
    if i < N - 1: s += ",\n "
  s+="]"
  if tostdout: print(s)
  return s


def check_transformation_matrices(geom, a, ir, jth,
    verbose=True, tol=1e-12) -> bool:
  """Transforms the metric to spherical KS and compares with analytic formula

  Test 1: covariant metric, gcov, at A = {ir, jth}
  1.1 sample gcov and Lambda_h2bl_cov at A
  1.2 transform gcov to gks using transofmration matrices
  1.3 compare to expected values at {r,th} at A

  Parameters
  ----------
  geom : dictionary
      nubhlight geom object
  a : Float
      dimensionless Kerr spin parameter
  ir   : Integer
      index of sample point in radial direction
  jth  : Integer
      index of sample point in angular theta-direction
  verbose : bool
      output steps to stdout
  tol : Float
      tolerance to relative error (wrt det g)

  Returns
  -------
  bool
      True if all checks passed

  Examples
  --------
    import hdf5_to_dict as io
    hdr = io.load_hdr("dump_00000010.h5")
    geom = io.load_geom(hdr,recalc=True)
    check_transformation_matrices(geom, -1, 64)
  """
  # sample gcov and h2bl at point A
  gcov_A = geom['gcov'][ir,jth]
  h2bl_A = geom['Lambda_h2bl_cov'][ir,jth]

  # sample r and theta, compute BL metric-related quantities
  r = geom['r'][ir,jth,0]; r2 = r*r
  a2 = a*a
  th= geom['th'][ir,jth,0]
  sth2= sin(th)**2
  Delta= r2 - 2*r + a2
  Sigma= r2 + a2*cos(th)**2
  A = (r2 + a2)**2 - a2*Delta*sin(th)**2

  if verbose:
    print ("r      = %19.11e" % r)
    print ("theta  = %19.11e" % th)
    print ("a      = %19.11e" % a)
    print ("Delta  = %19.11e" % Delta)
    print ("Sigma  = %19.11e" % Sigma)
    print ("A      = %19.11e" % A)

    # output metric
    print ("gcov_A = ")
    print_matrix (gcov_A)
    print ("")

    # output transformation matrix
    print ("h2bl_A = ")
    print_matrix (h2bl_A)
    print ("")

  # compute BL metric at A
  gks_A = np.zeros([4,4])

  for i in range(4):
   for j in range(4):
       for k in range(4):
        for l in range(4):
         gks_A[i,j] = gks_A[i,j] + h2bl_A[k,i]*h2bl_A[l,j]*gcov_A[k,l]
  if verbose:
    print ("gks_A = ")
    print_matrix (gks_A)
    print("")

  # expected values at {r, th}
  g_tt    = -1. + 2.*r/Sigma
  g_rr    =  1. + 2.*r/Sigma
  g_ff    = sth2*(Sigma + a2*g_rr*sth2)
  g_thth  = Sigma

  g_tr    = 2*r/Sigma
  g_tf    = -2*a*r*sth2/Sigma
  g_rf    = -a*g_rr*sth2
  det_g   = -Sigma**2*sth2

  if verbose:
    print ("Expected:")
    print (" g_tt    = %19.11e" % g_tt  )
    print (" g_rr    = %19.11e" % g_rr  )
    print (" g_thth  = %19.11e" % g_thth)
    print (" g_ff    = %19.11e" % g_ff  )
    print (" g_tr    = %19.11e" % g_tr  )
    print (" g_rf    = %19.11e" % g_rf  )
    print (" g_tf    = %19.11e" % g_tf  )
    print ("")

  # check gks_A
  gks_expected = np.array(
      [[ g_tt, g_tr, 0.0, g_tf],
       [ g_tr, g_rr, 0.0, g_rf],
       [ 0.0, 0.0, g_thth, 0.0],
       [ g_tf, g_rf, 0.0, g_ff]]
  )

  passed = True
  for i in range(4):
    for j in range(4):
      if abs(gks_A[i,j] - gks_expected[i,j])/abs(det_g) > tol:
        passed = False
        if verbose:
          print (f"WARNING: Significant mismatch in gks_A[{i},{j}]:")
          print ("  -- expected: %19.11e" % gks_expected[i,j])
          print ("  -- actual:   %19.11e" % gks_A[i,j])

  return passed
