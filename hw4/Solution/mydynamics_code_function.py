# 2022-4-16
# Created by LIU Ben for ME5114-homework4

import numpy as np
def axis2skew(w):
    if w.ndim == 1:
        return np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    else:
        (n,m)=w.shape
        if n == 1:
            w=w.transpose()
        return np.array([[0,-w[2,0],w[1,0]],[w[2,0],0,-w[0,0]],[-w[1,0],w[0,0],0]])

def exp_w(w,q):
    hatw=axis2skew(w)
    return np.identity(3)+hatw*np.sin(q)+hatw@hatw*(1-np.cos(q))

def exp_g(tw,q):
    w=tw[0:3]
    v=tw[3:6]
    if np.sum(w)==0:
        v=np.array([v]).transpose()
        return np.concatenate((np.concatenate((np.identity(3),v*q),axis=1),[[0,0,0,1]]))
    else:
        # finite pitch
        g=np.identity(4)
        g[0:3,0:3]=exp_w(w,q)
        g[0:3,3]=(np.identity(3)-exp_w(w,q))@(np.cross(w,v)+w*(w@v)*q)
        return g

def generate_twist(q1,q2,q3,w1,w2,w3):
    if w1+w2+w3==0:
        return np.array([w1,w2,w3,q1,q2,q3]).astype(float)
    else:
        return np.concatenate(([w1,w2,w3],np.cross(np.array([q1,q2,q3]),[w1,w2,w3]))).astype(float)

def generate_twist2(q,w):
    if np.sum(w)==0:
        # infinit pitch
        return np.concatenate((w,q))
    else:
        # finite pitch
        return np.concatenate((w,np.cross(q,w)))

def Ad_g(g):
    R=g[0:3,0:3]
    p=g[0:3,3]
    Ad=np.zeros((6,6))
    Ad[0:3,0:3]=R
    Ad[3:6,0:3]=axis2skew(p)@R
    Ad[3:6,3:6]=R
    return Ad

def Ad_g_inv(g):
    R=g[0:3,0:3]
    p=g[0:3,3]
    Ad=np.zeros((6,6))
    Ad[0:3,0:3]=R.T
    Ad[3:6,0:3]=-R.T@axis2skew(p)
    Ad[3:6,3:6]=R.T
    return Ad

def ad_V(tw):
    ad=np.zeros((6,6))
    w=tw[0:3]
    v=tw[3:6]
    ad[0:3,0:3]=axis2skew(w)
    ad[3:6,0:3]=axis2skew(v)
    ad[3:6,3:6]=axis2skew(w)
    return ad

def Jacobian_spatial(tw_all,q):
    # tw_all[:,i] is i-th twist
    J=tw_all.astype(float)
    if q.size == 1:
        return J
    else:
        for i in range(1,q.size):
            g=np.identity(4)
            for j in range(i):
                g=g@exp_g(tw_all[:,j],q[j])
            J[:,i]=Ad_g(g)@tw_all[:,i]
        return J
            
def Jacobian_body(tw_all,q,gst0):
    # tw_all[:,i] is i-th twist
    J=tw_all.astype(float)
    for i in range(q.size):
            g=np.identity(4)
            for j in range(i,q.size):
                g=g@exp_g(tw_all[:,j],q[j])
            J[:,i]=Ad_g_inv(g@gst0)@tw_all[:,i]
    return J


def xyz_spatial(x,y,z):
    return exp_w(np.array([0,0,1]),z)@exp_w(np.array([0,1,0]),y)@exp_w(np.array([1,0,0]),x)

def generate_SE3(Euler_zyx,p):
    # generate_SE3((theta1,theta2,theta3),p)
    R=xyz_spatial(Euler_zyx[0],Euler_zyx[1],Euler_zyx[2])
    T=np.identity(4)
    T[0:3,0:3]=R
    T[0:3,3]=p
    return T


def spatialinertia_fromCoM(T,m,I):
    # T = frame of oM frame expressed in spatial frame --- SE(3)
    # I is the inertia matrix written in CoM frame  --- (3 x 3)
    I_spatial=m*np.identity(6)
    I_spatial[0:3,0:3]=I
    I_spatial=Ad_g_inv(T).T@I_spatial@Ad_g_inv(T)
    return I_spatial
    
