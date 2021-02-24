import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from math import *
from scipy.signal import argrelextrema

import networkx as nx

import seaborn as sns
sns.set()
sns.set_context("notebook", font_scale=1.75, rc={"lines.linewidth": 4.0})

plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]



def readEdges(filename):
    edges=np.loadtxt(filename)
    edges.sort(axis=1)
    edges=edges[np.lexsort((edges[:, 1], edges[:, 0]))]
    return edges

def readTrigs(filename):
    trigs=np.loadtxt(filename)
    trigs.sort(axis=1)
    trigs=trigs[np.lexsort((trigs[:, 2], trigs[:, 1], trigs[:, 0]))]
    return trigs

def B1fromEdges(n, edges):
    B1=np.zeros((n, edges.shape[0]))

    for i in range(edges.shape[0]):
        B1[int(edges[i, 0]),i]=-1
        B1[int(edges[i, 1]),i]=1
    return B1

def B2fromTrig(n, edges, trigs):
    B2=np.zeros((edges.shape[0], trigs.shape[0]))

    for i in range(trigs.shape[0]):
        B2[np.where((edges==np.array([trigs[i, 0], trigs[i, 1]])).all(axis=1))[0][0] ,i]=1
        B2[np.where((edges==np.array([trigs[i, 0], trigs[i, 2]])).all(axis=1))[0][0] ,i]=-1
        B2[np.where((edges==np.array([trigs[i, 1], trigs[i, 2]])).all(axis=1))[0][0] ,i]=1
    return B2

def getRandomWeights(edges):
    return np.random.uniform(size=(edges.shape[0]))

def getAdjB1(B1):
    return np.diag(np.diag(B1.dot(B1.T)))-B1.dot(B1.T)

def getPosFromB1(B1):
    A=getAdjB1(B1)
    G = nx.from_numpy_matrix(np.array(A))  
    pos = nx.spring_layout(G)
    return np.array(list(pos.values()))

def myinv(W, thr=1e-12):
    w=np.diag(W)
    ans=np.zeros(w.shape)
    ans[np.abs(w)>thr]=1/w[np.abs(w)>thr]
    return np.diag(ans)

def getDt(B2, W, thr=1e-12):
    w=np.diag(W).reshape(-1, 1)
    on=np.ones((B2.shape[1], 1))
    tmp=np.multiply(on.dot(w.T), abs(B2.T))
    minval=np.zeros((tmp.shape[0]))
    for i in range(tmp.shape[0]):
       if np.sum(tmp[i, :]>thr) >= 3:
          minval[i]=np.min(tmp[i, tmp[i, :]>thr]) 
    Dt=np.diag(minval)
    return Dt

def getL1k(L1, thr=1e-8):
    return np.sum(np.abs(np.linalg.eigh(L1)[0])<thr)

def Sym(A):
    return 0.5*(A+A.T)

def scal(A,B):
    return np.trace(A.T @ B);

def getG_i(B1, B2, L1, x, w, e, eps, thr=1e-12):
    x = x.reshape(-1, 1)
    W = np.diag(np.sqrt(w)+eps*e)
    Dt = getDt(B2, W)
    P2 = myinv(W) @ B2 @ Dt @ Dt @ B2.T @ myinv(W)
    #P2 = myinv(W).dot(B2.dot(Dt.dot(Dt.dot(B2.T.dot(myinv(W))))))
    w_t = np.diag(W).reshape(-1, 1)
    on = np.ones((B2.shape[1], 1))
    tmp = np.multiply(on @ w_t.T, np.abs(B2.T));
    tmp[tmp==0]=2*np.sum(w_t)
    M = np.zeros(B2.T.shape);
    for i in range(M.shape[0]):
        M[i, np.argmin(tmp[i, :])]=1
    Add= np.diag(M.T @ np.diag(B2.T @ myinv(W) @ (x @ x.T) @ myinv(W) @ B2 @ Dt ))
    #Add = np.diag(M.T.dot(np.diag(B2.T.dot(myinv(W).dot(x.dot(x.T.dot(myinv(W).dot(B2.dot(Dt)))))))))
    Gi = 2.*Sym( B1.T @ B1 @ W @ (x @ x.T) - myinv(W) @ (x @ x.T) @ P2) + 2.*Add
    #Gi = 2.*Sym(B1.T.dot(B1.dot(W.dot(x.dot(x.T))))-myinv(W).dot(x.dot(x.T.dot(P2))))+2.*Add
    return 2.*(x.T @ L1 @ x)*Gi

def getDotE_con_new(B1, B2, w, e, eps, k, thr, alph, mu):
    matmask = np.diag(np.logical_not(np.abs(np.sqrt(w)+e*eps) < thr))
    E = np.diag(e)
    PE = np.multiply(E, matmask)

    L1_E = HodgeLW_fr(B1, B2, w, e, eps)
    vals, vecs=np.linalg.eigh(L1_E)
    idx = vals.argsort() 
    vals = vals[idx]
    vecs=vecs[:, idx]
    
    GE = np.zeros(L1_E.shape)
    for i in range(k+1):
        Gi=getG_i(B1, B2, L1_E, vecs[:,i], w, e, eps)
        GE = GE + Gi;
    
    L0=getL0(B1, w, e, eps)
    W=np.diag(np.sqrt(w)+eps*e)
    D12=np.diag(np.power(np.sum(np.abs(B1) @ W, 1), -1))
    vals, vecs=np.linalg.eigh(L0)
    idx = vals.argsort() 
    vals = vals[idx]
    x=vecs[:, 1].reshape(-1, 1)
    if vals[1]<mu:
        Add = 2.*B1.T @ D12 @ (x @ x.T) @ D12 @ B1 @ W - 2.*np.diag(np.diag(Sym(D12 @ (x @ x.T) @ L0)).reshape(-1, 1).T @ np.abs(B1))
        GE = GE - alph/mu*np.max([0, 1-vals[1]/mu])*Add
    
    
    kappa = scal(GE, PE)/scal(PE, PE)
    PGE = np.multiply(GE, matmask)
    dE = -PGE+kappa*PE
    return dE

def run_con_new(B1, B2, w, eps, e, beta, h, k, thr, alpha, mu, verbose, draw):
    e = e/np.linalg.norm(e)
    L1_E = HodgeLW_fr(B1, B2, w, e, eps)
    L0=getL0(B1, w, e, eps)
    track = [getFk_l2_con_new(L1_E, k, L0, alpha, mu)]
    t_cur = 0
    log = []
    ts = [0]
    h0=h
    tfin=1000
    nu=0.1
    i=0
    fl=0
    minval=1.0
    emin=np.zeros(e.shape)
    while (ts[-1]<tfin) and (i<=1000) :
        i+=1
        #print(i, end=' ')
        e0 = e
        fl=0
        while True and (fl==0):
            st=0
            while True:
                e = e0
                dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alpha, mu)
                E1 = np.diag(e)-h*dE
                e = np.diag(E1)
                e = e/np.linalg.norm(e)
                if verbose:
                    print('h: ', h, ' |  time: ', t_cur,'  ||   E_norm: ', round(scal(E1, E1), 3), ' , dE-orth: ', round(scal(dE, E1), 3), ' ||   F=', track[-1])
                if np.sum(np.sqrt(w)+eps*e < 0) > 0:
                    h = h/2
                    #print('!', h, end=' ')
                    #print("!", end=' ')
                    st+=1
                else:
                    break
                #if st>100:
                #    fl=1
                #    break
            L1_E = HodgeLW_fr(B1, B2, w, e, eps)
            L0=getL0(B1, w, e, eps)
            newval=getFk_l2_con_new(L1_E, k, L0, alpha, mu)
            if  newval < track[-1]:
                h = h*beta
                #print('W', h, end=' ')
                break
            else:
                h = h/beta  
                #print('X' ,h, end=' ')
            if h<1e-20 or (fl==1):
                h=h0
                check=0
                while check==0:
                    e = np.diag(E1)+nu/(i+1)*np.random.uniform(size=(e.shape[0]))
                    e = e/np.linalg.norm(e)
                    if np.sum(np.sqrt(w)+eps*e<0)>0:
                        check=1
                
                #print("jump, FUCK YOU!", ts[-1])
                break
        log.append(np.sort(np.linalg.eigvalsh(L1_E)))
        track.append(newval)
        if track[-1]<minval:
            minval=track[-1]
            emin=e
        #dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alpha, mu)
        t_cur =t_cur + h
        ts.append(t_cur)
    
    
    if draw:
        simpleDrawB1(B1, w, e, eps)
    

    return e, L1_E, ts, track, log, emin 


# In[19]:

def simpleDrawB1(B1, w, points, edges, eps=0, e=0):
    
    plt.figure(figsize=(8,8))
    plt.plot(points[:,0], points[:,1], 'o', color=colors[4], markersize=40)

    for i in range(points.shape[0]):
        plt.text(x=points[i,0], y=points[i,1], s=str(i), va='center', ha='center')

    for i in range(edges.shape[0]):
        plt.plot(points[edges[i].astype(int), 0], points[edges[i].astype(int), 1], color=colors[3], linewidth=12*((np.sqrt(w)+eps*e)**2)[i])
        plt.text(x=np.mean(points[edges[i].astype(int), 0]), 
                y=np.mean(points[edges[i].astype(int), 1]), s=str(round((np.sqrt(w)+eps*e)[i], 3)), va='center', ha='center',
                #rotation=180/np.pi*np.arctan(np.diff(points[edges[i].astype(int), 1])[0]/np.diff(points[edges[i].astype(int), 0])[0])
                color='grey', fontsize=14)
    plt.grid(False)
    plt.xticks([], [])
    plt.yticks([], [])
    plt.title(r'$\varepsilon=$'+str(round(eps, 3)))

def HodgeLW_fr(B1, B2, w, e=0, eps=0):
    W=np.diag(np.sqrt(w)+eps*e)
    Dt=getDt(B2, W)
    L1=W @ B1.T @ B1 @ W+myinv(W) @ B2 @ Dt @ Dt @ B2.T @ myinv(W)
    #L1=W.dot(B1.T.dot(B1.dot(W)))+myinv(W).dot(B2.dot(Dt.dot(Dt.dot(B2.T.dot(myinv(W))))))
    return L1

def getL0(B1, w, e, eps):
    W=np.diag(np.sqrt(w)+eps*e)
    D12=np.diag(np.power(np.sum(np.abs(B1) @ np.diag(w), 1), -1))
    L0=D12 @ B1 @ W @ W @ B1.T @ D12
    #L0=D12.dot(B1.dot(W.dot(W.dot(B1.T.dot(D12)))))
    return L0

def getFk_l2_con_new(L1, k, L0, alpha, mu, thr=1e-8):
    vals=np.sort(np.linalg.eigvalsh(L1))
    Fk=np.sum(np.power(vals[:k+1], 2))
    vals=np.sort(np.linalg.eigvalsh(L0))
    Fk+=.5*alpha*np.max([0, 1-vals[1]/mu])**2;
    return Fk
