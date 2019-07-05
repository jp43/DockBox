## this file corresponds to the file pyqcprot.pyx that can be cloned at https://github.com/synapticarbors/pyqcprot (or can be found in MDAnalysis package). Jordane Preto included function CenterCoords that is used to center coordinates before computing RMSD. No significant overhead with respect to the original version was reported.

# -----------------------------------------------------------------------------
#    Author(s) of Original Implementation:     
#                  Douglas L. Theobald
#                  Department of Biochemistry
#                  MS 009
#                  Brandeis University
#                  415 South St
#                  Waltham, MA  02453
#                  USA
#
#                  dtheobald@brandeis.edu
#                  
#                  Pu Liu
#                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
#                  665 Stockton Drive
#                  Exton, PA  19341
#                  USA
#
#                  pliu24@its.jnj.com
#
#                  For the original code written in C see:
#                  http://theobald.brandeis.edu/qcp/ 
#
#
#    Author of Python Port:
#                  Joshua L. Adelman
#                  Department of Biological Sciences
#                  University of Pittsburgh
#                  Pittsburgh, PA 15260
#                   
#                  jla65@pitt.edu
#
# 
#    If you use this QCP rotation calculation method in a publication, please
#    reference:
#
#      Douglas L. Theobald (2005)
#      "Rapid calculation of RMSD using a quaternion-based characteristic
#      polynomial."
#      Acta Crystallographica A 61(4):478-480.
#
#      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010)
#      "Fast determination of the optimal rotational matrix for macromolecular 
#      superpositions."
#      J. Comput. Chem. 31, 1561-1563. 
#
#
#    Copyright (c) 2009-2010, Pu Liu and Douglas L. Theobald 
#    Copyright (c) 2011       Joshua L. Adelman
#    All rights reserved.
#
#    Redistribution and use in source and binary forms, with or without modification, are permitted 
#    provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright notice, this list of 
#      conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list 
#      of conditions and the following disclaimer in the documentation and/or other materials 
#      provided with the distribution.
#    * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to 
#      endorse or promote products derived from this software without specific prior written 
#      permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------

import cython
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)

cdef double InnerProduct(np.ndarray[np.float64_t,ndim=1] A,
                 np.ndarray[np.float64_t,ndim=2] coords1,
                 np.ndarray[np.float64_t,ndim=2] coords2,
                 int N,
                 np.ndarray[np.float64_t,ndim=1] weight):
    """
    Calculate the inner product of two structures.
    InnerProduct(A, coords1, coords2, N, weight) --> (G1+G2)/2
    If weight array is not ``None``, calculate the weighted inner product.
            :Input:
                   - A[9]    -- inner product array (modified in place)
                   - coords1 -- reference structure
                   - coords2 -- candidate structure
                   - N       -- the size of the system
                   - weight  -- the weight array of size N: set to None if not needed
            :Output:
                   - A[9]    -- the inner product matrix
            :Returns:
                   - (G1 + G2) * 0.5; used as E0 in function :func:`FastCalcRMSDAndRotation`
            .. Warning::
                1. You MUST center the structures, coords1 and coords2, before calling this function.
                2. Please note how the structure coordinates are stored as 3xN arrays,
                   not Nx3 arrays as is also commonly used. The difference is
                   something like this for storage of a structure with 8 atoms::
                      Nx3: xyzxyzxyzxyzxyzxyzxyzxyz
                      3xN: xxxxxxxxyyyyyyyyzzzzzzzz
    """

    cdef double          x1, x2, y1, y2, z1, z2
    cdef unsigned int    i
    cdef double          G1, G2

    G1 = 0.0
    G2 = 0.0

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0

    if (weight != None):
        for i in xrange(N):
            x1 = weight[i] * coords1[0,i]
            y1 = weight[i] * coords1[1,i]
            z1 = weight[i] * coords1[2,i]

            G1 += x1*coords1[0,i] + y1*coords1[1,i] + z1*coords1[2,i]

            x2 = coords2[0,i]
            y2 = coords2[1,i]
            z2 = coords2[2,i]

            G2 += weight[i] * (x2*x2 + y2*y2 + z2*z2)

            A[0] +=  (x1 * x2)
            A[1] +=  (x1 * y2)
            A[2] +=  (x1 * z2)

            A[3] +=  (y1 * x2)
            A[4] +=  (y1 * y2)
            A[5] +=  (y1 * z2)

            A[6] +=  (z1 * x2)
            A[7] +=  (z1 * y2)
            A[8] +=  (z1 * z2)

    else:
        for i in xrange(N):
            x1 = coords1[0,i]
            y1 = coords1[1,i]
            z1 = coords1[2,i]

            G1 += (x1*x1 + y1*y1 + z1*z1)

            x2 = coords2[0,i]
            y2 = coords2[1,i]
            z2 = coords2[2,i]

            G2 += (x2*x2 + y2*y2 + z2*z2)

            A[0] +=  (x1 * x2)
            A[1] +=  (x1 * y2)
            A[2] +=  (x1 * z2)

            A[3] +=  (y1 * x2)
            A[4] +=  (y1 * y2)
            A[5] +=  (y1 * z2)

            A[6] +=  (z1 * x2)
            A[7] +=  (z1 * y2)
            A[8] +=  (z1 * z2)

    return (G1 + G2) * 0.5

cdef double FastCalcRMSDAndRotation(np.ndarray[np.float64_t,ndim=1] rot, np.ndarray[np.float64_t,ndim=1] A, double E0, int N):
    """
    Calculate the RMSD, and/or the optimal rotation matrix.
    FastCalcRMSDAndRotation(rot, A, E0, N)
            :Input:
                    - rot[9]  -- rotation matrix (modified in place)
                    - A[9]    -- the inner product of two structures
                    - E0      -- (G1 + G2) * 0.5
                    - N       -- the size of the system
            :Output:
                    - rot[9]   -- the rotation matrix in the order of xx, xy, xz, yx, yy, yz, zx, zy, zz
                    - rmsd     -- the RMSD value
            :Returns:
                    - only the rmsd was calculated if rot is None
                    - both the RMSD & rotational matrix calculated if rot is not None
    """
    cdef double rmsd
    cdef double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
    cdef double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
    cdef double SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
    cdef double SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
    cdef double SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy

    cdef np.ndarray[np.float64_t,ndim=1] C = np.zeros(4, dtype=A.dtype)
    cdef unsigned int i
    cdef double mxEigenV
    cdef double oldg = 0.0
    cdef double b, a, delta, rms, qsqr
    cdef double q1, q2, q3, q4, normq
    cdef double a11, a12, a13, a14, a21, a22, a23, a24
    cdef double a31, a64, a33, a34, a41, a42, a43, a44
    cdef double a2, x2, y2, z2
    cdef double xy, az, zx, ay, yz, ax
    cdef double a3344_4334, a6444_4234, a6443_4233, a3143_4133,a3144_4134, a3142_4164
    cdef double evecprec = 1e-6
    cdef double evalprec = 1e-14

    cdef double a1644_1423, a1224_1422, a1223_1642, a1124_1421, a1123_1641, a1122_1221
    Sxx = A[0]
    Sxy = A[1]
    Sxz = A[2]
    Syx = A[3]
    Syy = A[4]
    Syz = A[5]
    Szx = A[6]
    Szy = A[7]
    Szz = A[8]

    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    C[0] = (Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz)))

    mxEigenV = E0
    for i in xrange(50):
        oldg = mxEigenV
        x2 = mxEigenV*mxEigenV
        b = (x2 + C[2])*mxEigenV
        a = b + C[1]
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a))
        mxEigenV -= delta
        if (fabs(mxEigenV - oldg) < fabs((evalprec)*mxEigenV)):
            break

    #if (i == 50):
    #   print "\nMore than %d iterations needed!\n" % (i)

    # the fabs() is to guard against extremely small, but *negative* numbers due to npfloat point error
    rms = sqrt(fabs(2.0 * (E0 - mxEigenV)/N))

    if (rot == None):
        return rms # Don't bother with rotation.

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx
    a31 = a13; a64 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV
    a3344_4334 = a33 * a44 - a43 * a34; a6444_4234 = a64 * a44-a42*a34
    a6443_4233 = a64 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4164 = a31 * a42-a41*a64
    q1 =  a22*a3344_4334-a23*a6444_4234+a24*a6443_4233
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
    q3 =  a21*a6444_4234-a22*a3144_4134+a24*a3142_4164
    q4 = -a21*a6443_4233+a22*a3143_4133-a23*a3142_4164

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

# The following code tries to calculate another column in the adjoint matrix when the norm of the
#   current column is too small.
#   Usually this commented block will never be activated.  To be absolutely safe this should be
#   uncommented, but it is most likely unnecessary.

    if (qsqr < evecprec):
        q1 =  a12*a3344_4334 - a13*a6444_4234 + a14*a6443_4233
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
        q3 =  a11*a6444_4234 - a12*a3144_4134 + a14*a3142_4164
        q4 = -a11*a6443_4233 + a12*a3143_4133 - a13*a3142_4164
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

        if (qsqr < evecprec):
            a1644_1423 = a13 * a24 - a14 * a23
            a1224_1422 = a12 * a24 - a14 * a22
            a1223_1642 = a12 * a23 - a13 * a22
            a1124_1421 = a11 * a24 - a14 * a21
            a1123_1641 = a11 * a23 - a13 * a21
            a1122_1221 = a11 * a22 - a12 * a21

            q1 =  a42 * a1644_1423 - a43 * a1224_1422 + a44 * a1223_1642
            q2 = -a41 * a1644_1423 + a43 * a1124_1421 - a44 * a1123_1641
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
            q4 = -a41 * a1223_1642 + a42 * a1123_1641 - a43 * a1122_1221
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

            if (qsqr < evecprec):
                q1 =  a64 * a1644_1423 - a33 * a1224_1422 + a34 * a1223_1642
                q2 = -a31 * a1644_1423 + a33 * a1124_1421 - a34 * a1123_1641
                q3 =  a31 * a1224_1422 - a64 * a1124_1421 + a34 * a1122_1221
                q4 = -a31 * a1223_1642 + a64 * a1123_1641 - a33 * a1122_1221
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4

                if (qsqr < evecprec):
                    # if qsqr is still too small, return the identity matrix. #
                    rot[0] = rot[4] = rot[8] = 1.0
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0

                    return rms


    normq = sqrt(qsqr)
    q1 /= normq
    q2 /= normq
    q3 /= normq
    q4 /= normq

    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4

    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2

    rot[0] = a2 + x2 - y2 - z2
    rot[1] = 2 * (xy + az)
    rot[2] = 2 * (zx - ay)
    rot[3] = 2 * (xy - az)
    rot[4] = a2 - x2 + y2 - z2
    rot[5] = 2 * (yz + ax)
    rot[6] = 2 * (zx + ay)
    rot[7] = 2 * (yz - ax)
    rot[8] = a2 - x2 - y2 + z2

    return rms

cdef void CenterCoords(np.ndarray[np.float64_t,ndim=2] coords, int N, np.ndarray[np.float64_t,ndim=1] weights):

    cdef double          xsum, ysum, zsum, wsum
    cdef unsigned int    i

    xsum = ysum = zsum = 0.0

    if (weights != None):
        wsum = 0.0
        for i in xrange(N):
            xsum += weights[i] * coords[0,i]
            ysum += weights[i] * coords[1,i]
            zsum += weights[i] * coords[2,i]

            wsum += weights[i]

        xsum /= wsum
        ysum /= wsum
        zsum /= wsum

    else:
        for i in xrange(N):
            xsum += coords[0,i]
            ysum += coords[1,i]
            zsum += coords[2,i]

        xsum /= N
        ysum /= N
        zsum /= N

    for i in xrange(N):
        coords[0,i] -= xsum
        coords[1,i] -= ysum 
        coords[2,i] -= zsum

    return

@cython.boundscheck(False)
@cython.wraparound(False)
def CalcRMSDRotationalMatrix(np.ndarray[np.float64_t,ndim=2] ref,
                             np.ndarray[np.float64_t,ndim=2] conf,
                             np.ndarray[np.float64_t,ndim=1] rot,
                             np.ndarray[np.float64_t,ndim=1] weights):
    """
    Calculate the RMSD & rotational matrix.
    CalcRMSDRotationalMatrix(ref, conf, N, rot, weights):
            :Input:
                   - ref     -- reference structure coordinates (*must* be `numpy.float64`)
                   - conf    -- candidate structure coordinates (*must* be `numpy.float64`)
                   - rot[9]  -- array to store rotation matrix; set to None if only calculating rmsd (modified in place)
                   - weight  -- the weight array of size len; set to None if not needed
            :Output:
                   - rot[9]  -- rotation matrix
            :Returns:
                   - RMSD value
    .. Note:: All arrays *must* be of type `numpy.float64`.
    """
    cdef double rmsd
    cdef int N = conf.shape[1]
    cdef double E0
    cdef np.ndarray[np.float64_t,ndim=1] A = np.zeros(9, dtype=ref.dtype)
 
    CenterCoords(ref, N, weights)
    CenterCoords(conf, N, weights)

    E0 = InnerProduct(A, conf, ref, N, weights)
    rmsd = FastCalcRMSDAndRotation(rot, A, E0, N)

    return rmsd

    # the following lines are here to show how to apply the rotational matrix in order to recover the right RMSD
    #cdef double trmsd = 0.0
    #cdef np.ndarray[np.float64_t,ndim=1] trot = np.zeros(3, dtype=ref.dtype)

    #rmsd = CalcRMSDRotationalMatrix(ref, conf, rot, weights)

    #for k in xrange(N):
    #    for i in xrange(3):
    #        trot[i] = 0.0
    #        for j in xrange(3):
    #            trot[i] += rot[3*i+j] * ref[j,k]
    #    trmsd += (trot[0]-conf[0,k])**2 + (trot[1]-conf[1,k])**2 + (trot[2]-conf[2,k])**2
    #trmsd = sqrt(trmsd/N)
    #print rmsd, trmsd
