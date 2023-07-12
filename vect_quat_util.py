#!/usr/bin/python
import numpy as np
import math as m
from numba import jit, float64
from subtool import *

def vector_dot(u,v):
    return np.dot(u,v)
def vector_norm(u):
    return np.linalg.norm(u)
def vector_angle(u,v):
    t = vector_dot(u,v) / (vector_norm(u)*vector_norm(v))
    if (t > 1):
        t = 1
    elif (t < -1):
        t = -1
    return m.acos(t)
    #return m.acos(vector_dot(u,v) / (vector_norm(u)*vector_norm(v)))
def vector_cross(u,v):
    s = np.zeros((3,))
    s[0] = u[1]*v[2] - u[2]*v[1]
    s[1] = u[2]*v[0] - u[0]*v[2]
    s[2] = u[0]*v[1] - u[1]*v[0]
    return s

def quat_normalize(q):
    n = m.sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
    return [q[0], q[1]/n, q[2]/n, q[3]/n]

def quat_multiply(p,q): # out = p*q
    dest = np.zeros(4)
    dest[0] =  p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3];
    dest[1] =  p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2];
    dest[2] =  p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1];
    dest[3] =  p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0];
    return dest

def quat_norm(q): # return |q|
    return m.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])

def quat_dot(p,q): # p dot q
  return p[0]*q[0]+p[1]*q[1]+p[2]*q[2]+p[3]*q[3];

def quat_conj(q): # return conj(q)
    return [q[0], -q[1], -q[2], -q[3]]

def quat_vec_rot(src,q,dest=None):
  # quat_vec_rot(vector &dest, vector &src, quaternion &q) dest = q*src*conj(q)
  if dest is not None:
    return _quat_vec_rot(src,q,dest)
  dest = np.zeros(3)
  aa=[q[0]*q[0], q[1]*q[1], q[2]*q[2], q[3]*q[3]]
  ab=[q[0]*q[1], q[0]*q[2], q[0]*q[3], q[1]*q[2], q[1]*q[3], q[2]*q[3]]
  dest[0] = (aa[0]+aa[1]-aa[2]-aa[3])*src[0]+\
            ((ab[3]-ab[2])*src[1]+(ab[1]+ab[4])*src[2])*2.0;
  dest[1] = (aa[0]-aa[1]+aa[2]-aa[3])*src[1]+\
            ((ab[2]+ab[3])*src[0]+(ab[5]-ab[0])*src[2])*2.0;
  dest[2] = (aa[0]-aa[1]-aa[2]+aa[3])*src[2]+\
            ((ab[4]-ab[1])*src[0]+(ab[0]+ab[5])*src[1])*2.0;
  return dest;

#@jit(float64(float64,float64,float64,float64,float64),nopython=True,nogil=True,parallel=True)
#@jit
def _quat_vec_rot(src,q,dest):
  """Vectorized version of quat_vec_rot.
  tmp1 should be passed as np.zeros(4), tmp2 should be passed as np.ones(6)
  """
  tmp1 = np.zeros(4)
  tmp2 = np.ones(6)
  tmp1 = q*q
  
  tmp2[:3] *= q[0]
  tmp2[[0,3,4]] *= q[1]
  tmp2[[1,3,5]] *= q[2]
  tmp2[[2,4,5]] *= q[3]
  dest[0] = (tmp1[0]+tmp1[1]-tmp1[2]-tmp1[3])*src[0]+\
            ((tmp2[3]-tmp2[2])*src[1]+(tmp2[1]+tmp2[4])*src[2])*2.0;
  dest[1] = (tmp1[0]-tmp1[1]+tmp1[2]-tmp1[3])*src[1]+\
            ((tmp2[2]+tmp2[3])*src[0]+(tmp2[5]-tmp2[0])*src[2])*2.0;
  dest[2] = (tmp1[0]-tmp1[1]-tmp1[2]+tmp1[3])*src[2]+\
            ((tmp2[4]-tmp2[1])*src[0]+(tmp2[0]+tmp2[5])*src[1])*2.0;

  return dest

def quat_fvu_rot(fvu,q):
    f = quat_vec_rot(fvu[0],q)
    v = quat_vec_rot(fvu[1],q)
    u = quat_vec_rot(fvu[2],q)
    return np.array([f,v,u])

def tu2rotquat(theta,u,dest=None): # theta and u to rotational quat
    u = np.divide(u,vector_norm(u))
    if dest is None:
      dest = np.zeros(4)
    dest[0] = m.cos(theta/2.0)
    dest[1] = u[0]*m.sin(theta/2.0)
    dest[2] = u[1]*m.sin(theta/2.0)
    dest[3] = u[2]*m.sin(theta/2.0)

    return dest

def rotquat2tu(q): # rotational quat to theta and u
    theta = 2.0*m.acos(q[0])
    u = [0]*3
    u[0] = q[1] / m.sin(theta/2.0)
    u[1] = q[2] / m.sin(theta/2.0)
    u[2] = q[3] / m.sin(theta/2.0)
    return theta,u


