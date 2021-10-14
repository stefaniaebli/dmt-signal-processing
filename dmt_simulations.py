#!/usr/bin/env python3
import numpy as np
import random
import gudhi as gd
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import dmtenergy as dmt
import dmtvisual as dmtvis



X=np.load('./figures/X.npy',allow_pickle=True)
kX=np.load('./figures/kX.npy',allow_pickle=True)
s1=np.load('./figures/signal.npy')

total_l=[]
total_l_rand=[]
for j in range(5):

    steps=[5,10,20,30,50,70,100]
    optimal_simul=dmt.simulation_collapses(X=X,kX=kX,dimq=1,signal=s1,steps=steps,random=False)
    random_simul=dmt.simulation_collapses(X=X,kX=kX,dimq=1,signal=s1,steps=steps,random=True)
    total_l.append(optimal_simul)
    total_l_rand.append(random_simul)
    print(j)
    
np.save('./figures/loos_optimal.npy',total_l)
np.save('./figures/loos_random.npy',total_l_rand)