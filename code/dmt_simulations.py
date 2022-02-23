#!/usr/bin/env python3
import numpy as np
import random
import gudhi as gd
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import dmtsignal as dmt
import dmtvisual as dmtvis
import numpy as np


X=np.load('../data/X.npy',allow_pickle=True)
kX=np.load('../data/kX.npy',allow_pickle=True)
s1=np.load('../data/signal.npy')
#s1=np.load("../data/random_signal.npy",allow_pickle=True)
#s1=np.load("../data/normal_signal.npy",allow_pickle=True)
#s1=np.load("../data/center_signal.npy",allow_pickle=True)

lenc=len(X[2])
for k in range(1):
    #s1=np.random.random(len(X[1]))
    total_l=[]
    total_l_rand=[]
    for j in range(10):
        
        all_X,collapses,all_losses,total_loss,all_signals,phispsis= dmt.sequence_optimal_up_collapses(X=X,kX=kX,dimq=1,signal=s1,steps=lenc,random=False)
        all_X_rand,collapses_rand,all_losses_rand,total_loss_rand,all_signals_rand,phispsis_rand= dmt.sequence_optimal_up_collapses(X=X,kX=kX,dimq=1,signal=s1,steps=lenc,random=True)
 
        total_l.append([total_loss,all_losses])
        total_l_rand.append([total_loss_rand,all_losses_rand])
        print(j)
    print ("step:", k)  
    #np.save('../figures/data_optimal_sim{}.npy'.format(k),total_l)
    #np.save('../figures/data_random_sim{}.npy'.format(k),total_l_rand)
    np.save('../data/data_optimal_height_sim{}.npy'.format(k),total_l)
    np.save('../data/data_random_height_sim{}.npy'.format(k),total_l_rand)