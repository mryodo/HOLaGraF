import seaborn as sns
import networkx as nx
from scipy.signal import argrelextrema
from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import cm

from time import time

sns.set()
sns.set_context("notebook", font_scale=1.75, rc={"lines.linewidth": 4.0})

plt.style.use(
    'https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

from copy import deepcopy


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


def getRandomWeights(edges):
    return np.random.uniform(size=(edges.shape[0]))


class Graph:
    colors= plt.rcParams["axes.prop_cycle"].by_key()["color"]

    def __init__(self, n, edges, trigs, w, eps0, e):
        self.n, self.edges, self.trigs, self.w, self.eps0, self.e=n, edges, trigs, w, eps0, e
        self.B1=self.B1fromEdges(n, edges)
        self.B2=self.B2fromTrig(n, edges, trigs)

    def B1fromEdges(self, n, edges):
        B1=np.zeros((n, edges.shape[0]))

        for i in range(edges.shape[0]):
            B1[int(edges[i, 0]),i]=-1
            B1[int(edges[i, 1]),i]=1
        return B1

    def B2fromTrig(self, n, edges, trigs):
        B2=np.zeros((edges.shape[0], trigs.shape[0]))

        for i in range(trigs.shape[0]):
            B2[np.where((edges==np.array([trigs[i, 0], trigs[i, 1]])).all(axis=1))[0][0] ,i]=1
            B2[np.where((edges==np.array([trigs[i, 0], trigs[i, 2]])).all(axis=1))[0][0] ,i]=-1
            B2[np.where((edges==np.array([trigs[i, 1], trigs[i, 2]])).all(axis=1))[0][0] ,i]=1
        return B2

    def getW(self):
        return np.diag(np.sqrt(self.w)+self.eps0*self.e)

    def getAdjB1(self):
        return np.diag(np.diag(self.B1 @ self.B1.T))-self.B1 @ self.B1.T
    
    def getAdjB1W(self):
        W=self.getW()
        return np.diag(np.diag(self.B1 @ W @ W @ self.B1.T))-self.B1 @ W @ W @ self.B1.T

    def getPositions(self):
        A=self.getAdjB1()
        temp = nx.from_numpy_matrix(np.array(A))  
        pos = nx.spring_layout(temp)
        return np.array(list(pos.values()))

    def simpleDrawB1(self):
        points=self.getPositions()
        realW=np.diag(self.getW()**2)
        plt.figure(figsize=(8,8))
        plt.plot(points[:,0], points[:,1], 'o', color=self.colors[4], markersize=40)

        for i in range(points.shape[0]):
            plt.text(x=points[i,0], y=points[i,1], s=str(i), va='center', ha='center')

        for i in range(self.edges.shape[0]):
            plt.plot(points[self.edges[i].astype(int), 0], points[self.edges[i].astype(int), 1], color=self.colors[3], linewidth=12*realW[i])
            plt.text(x=np.mean(points[self.edges[i].astype(int), 0]), 
                    y=np.mean(points[self.edges[i].astype(int), 1]), s=str(round(realW[i], 3)), va='center', ha='center',
                    #rotation=180/np.pi*np.arctan(np.diff(points[edges[i].astype(int), 1])[0]/np.diff(points[edges[i].astype(int), 0])[0])
                    color='grey', fontsize=14)
        plt.grid(False)
        plt.xticks([], [])
        plt.yticks([], [])
        plt.title(r'$\varepsilon=$'+str(round(self.eps0, 3)))


def inv(A, thr=1e-10):
    w=np.diag(A)
    ans=np.zeros(w.shape[0])
    ans=1./w
    return np.diag(ans)

def Sym(A):
    return 0.5*(A+A.T)

def scal(A, B):
    return np.trace(A.T @ B)

def getDt(G, W):
    w_t = np.diag(W).reshape(-1, 1)
    ones2 = np.ones(G.B2.shape[1]).reshape(-1, 1)
    tmp = np.multiply(ones2 @ w_t.T, np.abs(G.B2).T)
    minval=np.zeros(tmp.shape[0])
    for i in range(tmp.shape[0]):
        if (tmp[i, np.nonzero(tmp[i])]).shape[1]==3:
            minval[i]=np.min(tmp[i, np.nonzero(tmp[i])])
    return np.diag(minval)

def getL1(G: Graph):
    W = G.getW()
    Dt = getDt(G, W)
    L1 = W @ G.B1.T @ G.B1 @ W + inv(W) @ G.B2 @ Dt @ Dt @ G.B2.T @ inv(W)
    return L1

def getL0(G: Graph):
    W=G.getW()
    A=G.getAdjB1W()+np.eye(G.B1.shape[0])
    d= (A @ np.ones((A.shape[0], 1))).flatten()
    L0=np.diag(d)-A
    L0=np.diag(np.power(d, -1/2).flatten()) @ L0 @ np.diag(np.power(d, -1/2).flatten())
    return L0

def getFk_l2(G, k, p , thrs):
    return getFk1(G, k, p , thrs)+getFk2(G, k, p , thrs)

def getFk1(G, k, p , thrs):
    L1=getL1(G)
    vals, _=np.linalg.eigh(L1)
    vals.sort()
    return 0.5*np.sum(np.power(vals[:k+p+1], 2))

def getFk2(G, k, p , thrs):
    L0=getL0(G)
    vals, _=np.linalg.eigh(L0)
    vals.sort()
    return 0.5*thrs['alpha']*np.max([ 0, 1-vals[1]/thrs['mu'] ])**2


def getNumGrad(G: Graph, k, p, thrs, thr0, normcor):
    gradNum=np.zeros(G.w.shape[0])
    delt=1e-4
    
    initValue=getFk_l2(G, k, p, thrs)
    W=G.getW()
    for i in range(G.w.shape[0]):
        if np.diag(np.abs(W))[i]>thr0:
            move=np.zeros(G.w.shape[0])
            move[i]=delt
            G1=deepcopy(G)
            G1.e=G1.e+move
            newValue=getFk_l2(G1, k, p, thrs)
            gradNum[i]=(newValue-initValue)/delt

    gradNum=gradNum/G.eps0

    if normcor:
        mask=(np.diag(W)>thr0)
        PE=np.multiply(mask, e)
        GPE=np.multiply(mask, gradNum)
        kappa=-(GPE.T @ PE) / (PE.T @ PE)
        gradNum=GPE+kappa*PE

    return gradNum


def getEigSort(A):
    vals, vecs=np.linalg.eigh(A)
    idx = np.argsort(vals)
    vals=vals[idx]
    vecs=vecs[:, idx]
    return vals, vecs

def getM(G: Graph, W, Dt):
    M=np.zeros(G.B2.T.shape)
    for i in range(M.shape[0]):
        ind=np.where( np.abs( np.diag(W) - Dt[i,i])<1e-08)[0][0]
        M[i, ind]=1
    return M


def getGrad(G: Graph, k: int, p: int, thrs: dict, thr0, normcor):
    L1=getL1(G)
    vals, vecs=getEigSort(L1)

    W=G.getW()
    Dt=getDt(G, W)
    P2=inv(W) @ G.B2 @ Dt @ Dt @ G.B2.T @ inv(W)
    preMat=G.B1.T @ G.B1 @ W
    preMat2=Dt @ G.B2.T @ inv(W)
    preMat3=inv(W) @ G.B2
    M=getM(G, W, Dt)

    grad=np.zeros(G.e.shape)
    for i in range(k+p+1):
        x=vecs[:, i].reshape(-1, 1)
        lam=vals[i]
        dE1=2*lam*Sym(preMat @ (x @ x.T) - P2 @ (x @ x.T) @ inv(W))
        A=preMat2 @ (x @ x.T) @ preMat3 
        dE2=2*lam*np.diag(M.T @ np.diag(A))
        grad+=np.diag(dE1+dE2)

    L0=getL0(G)
    vals, vecs=getEigSort(L0)
    lam=vals[1]
    x=vecs[:, 1].reshape(-1, 1)

    A=G.getAdjB1W()+np.eye(G.B1.shape[0])
    d= (A @ np.ones((A.shape[0], 1))).flatten()

    if 1- lam/thrs["mu"]>0:
        C=-thrs["alpha"]/thrs["mu"]*np.max([0, 1- lam/thrs["mu"]])
        #dE3=2*np.diag(W @ G.B1.T @ np.diag(Sym(np.diag(1./d) @ (x @ x.T) @ L0))) @ np.diag(G.B1.T @ np.ones((G.B1.T.shape[1] ,1)))
        tmp=np.multiply(W @ G.B1.T, G.B1.T) @ np.diag(Sym(np.diag(np.power(d, -1.0)) @ (x @ x.T) @ L0 )).reshape(-1, 1)  
        dE4=-2*np.diag(tmp.flatten())
        dE5=2*G.B1.T @ np.diag(np.power(d, -1/2)) @ (x @ x.T) @ np.diag(np.power(d, -1/2)) @ G.B1 @ W

        #dE6=-2*np.diag(W @ G.B1.T @ W @ np.diag(np.diag(np.power(d, -1/2)) @ (x @ x.T) @ np.diag(np.power(d, -1/2)))) @ np.diag(G.B1.T @ np.ones((G.B1.T.shape[1] ,1)))

        #grad += np.diag(dE3+dE4+dE5+dE6)
        grad += C*np.diag(dE4+dE5)

    
    if normcor:
        mask=(np.diag(W)>thr0)
        PE=np.multiply(mask, e)
        GPE=np.multiply(mask, grad)
        kappa=-(GPE.T @ PE) / (PE.T @ PE)
        grad=GPE+kappa*PE


    return grad


edges=readEdges('8.edges')
trigs=readTrigs('8.trigs')
n=8
#w=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
w=np.array([1.55, 0.86, 0.22, 1.62, 0.3, 1.02, 0.69, 1.87, 0.71, 1.57])
w=getRandomWeights(edges)+0.1


eps0=0.1
e0=-np.ones(w.shape)
e0=-np.random.uniform(size=w.shape)
e0/=np.linalg.norm(e0)
#e0=-np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0])


print(np.sqrt(w)+eps0*e0)
G=Graph(n, edges, trigs, w, eps0, e0)


G.simpleDrawB1()
plt.savefig('test.png', format='png')

thrs={
    'alpha': 10,
    'mu': 1.0
}

k=1
p=0
thr0=0


st_time=time()
numGrad=getNumGrad(G, k, p, thrs, thr0, 0)
fin_time1=time()
Grad=getGrad(G, k, p, thrs, thr0, 0)
fin_time2=time()
print(numGrad-Grad)
print("Time numerical:", fin_time1-st_time)
print("Time theoretical:", fin_time2-fin_time1)


