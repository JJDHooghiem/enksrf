# Example using kronecker matmul 
# Copyright (C) 2024 J.J.D. Hooghiem 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np 
import pyenkf 
import time as tt

# Example of kron(A1,A2) @ B
# usefull when the covariance matrix can be split into
# kronecker factors

# setup random matrices
# note that for large N and M the numpy equivalent 
# will try to allocate a lot of memory since the 
# full matrix is computed
#
# current implementation focussus on M of the order of a couple of 100
# and N is a couple of 1000
N=10
space=np.zeros((N,N),dtype=np.float64,order='F')
M=400
time=np.zeros((M,M),dtype=np.float64,order='F')
a=np.random.randn(N,N)
space[:,:]=a @ a.T

a=np.random.randn(M,M)
time[:,:] =a @ a.T
#
n=150
rn=np.random.randn(len(time)*len(space),n)

# compute the cholesky decomposition
# gives two lower triangular matrices
ti=np.linalg.cholesky(time)
sp=np.linalg.cholesky(space)

# dtrkmm as distributed with pyenkf
res=np.zeros((M*N,n),dtype=np.float64,order='F')
t=tt.time()
pyenkf.enkf_core.dtrkmm(sp,ti,rn,res)
print(tt.time()-t)


## reference routine 
if N*M<1000:
    t=tt.time()
    tot=np.kron(space,time)
    CD=np.linalg.cholesky(tot)
    cor=np.dot(CD,rn)
    print(tt.time()-t)
    print(np.allclose(res,cor))

# most efficient numpy
# skip if is to large (memory and speed)
if N*M<=10000:
    t=tt.time()
    cor=np.dot(np.kron(sp,ti),rn)
    print(tt.time()-t)
    print(np.allclose(res,cor))

exit()
