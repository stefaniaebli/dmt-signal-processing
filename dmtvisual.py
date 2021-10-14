#!/usr/bin/env python3
import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix
import gudhi as gd
import copy
import random
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.cm as cm
from dmtenergy import *


###Plotting Function
def get_positions(points,simplices, dim):
    polygons = list()
    for i, simplex in enumerate(simplices[dim].keys()):
        assert simplices[dim][simplex] == i  # Dictionary is ordered.
        polygon = list()
        for vertex in simplex:
            polygon.append(points[vertex])
        polygons.append(polygon)
    return polygons
      
def value2color(values):
    values -= values.min()
    values = values/(values.max()-values.min())
    return mpl.cm.viridis(values)

def plot_nodes(colors,points, ax=None,**kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    ax.scatter(points[:, 0], points[:, 1], c=colors, **kwargs)
    return ax.figure, ax

def plot_edges(colors,points,simplices, ax=None, **kwargs):   
    lines = get_positions(points,simplices, 1)
    if ax is None:
        fig, ax = plt.subfigs()
    colors = value2color(colors)
    collection = mpl.collections.LineCollection(lines, colors=colors, **kwargs)
    ax.add_collection(collection)
    ax.autoscale()

def plot_edges_normed(colors,points,simplices, ax=None, **kwargs):   
    lines = get_positions(points,simplices, 1)
    if ax is None:
        fig, ax = plt.subfigs()
    colors = mpl.cm.viridis(colors)
    collection = mpl.collections.LineCollection(lines, colors=colors, **kwargs)
    ax.add_collection(collection)
    ax.autoscale()

def plot_triangles(colors,points,simplices, ax=None, **kwargs):
    triangles = get_positions(points,simplices, 2)
    if ax is None:
        fig, ax = plt.subfigs()
    colors = value2color(colors)
    for triangle, color in zip(triangles, colors):
        triangle = plt.Polygon(triangle, color=color, **kwargs)
        ax.add_patch(triangle)
    ax.autoscale()
    
def plot_triangles_plain(color,points, simplices, ax=None, **kwargs):
    triangles = get_positions(points,simplices, 2)
    if ax is None:
        fig, ax = plt.subfigs()
    if type(color)==str:   
        colors = len(simplices[1])*[color]
    if type(color)!=str:   
        colors =color
    
    for triangle, color in zip(triangles, colors):
        triangle = plt.Polygon(triangle, color=color, **kwargs)
        ax.add_patch(triangle)
    ax.autoscale()
    
def plot_edges_plain(colors,points,simplices, ax=None, **kwargs):   
    lines = get_positions(points,simplices, 1)
    if ax is None:
        fig, ax = plt.subfigs()
    collection = mpl.collections.LineCollection(lines, colors=colors, **kwargs)
    ax.add_collection(collection)
    ax.autoscale()


    
def from_indices_to_color(X,collapses,dimq):
    down=[]
    up=[]
    for collapse in collapses:
        down.append(X[dimq][frozenset(collapse[0])])
        up.append(X[dimq+1][frozenset(collapse[1])])

    return(down,up)    
    

    
def plot_sequence_edgecollapses(X,dimq,signal,collapsed_X,collapsed_signal,phipsis,collapses,points,color_tri,size_nodes,size_lines,type_collapse="up",tot_max=None,tot_min=None):
    s0 = ['black']*len(X[0])#np.zeros(len(simplices[0]))
    
    s1 = signal
    s1phi=phipsis
    s1loss=np.abs(s1-s1phi)
    s1c=collapsed_signal
    a=[s1,s1phi,s1loss,s1c]
    if tot_min==tot_max==None:
       
        tot_max=max(map(max, a))
        tot_min=min(map(min, a))
    
    
    s1 = (s1-tot_min)/(tot_max-tot_min)
    s1phi=(s1phi-tot_min)/(tot_max-tot_min)
    s1loss=(s1loss-tot_min)/(tot_max-tot_min)
    s1c=(s1c-tot_min)/(tot_max-tot_min)
    
    
   
    ticks=np.around(np.arange(tot_min,tot_max,0.2),1)
    
    
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20, 20))
    plot_nodes(s0, points,ax=axes[0,0], zorder=3,s=size_nodes)
    plot_edges_normed(s1,points,X, ax=axes[0,0], zorder=2,linewidths=size_lines)
    if len(X)>2: plot_triangles_plain(color_tri,points,X, ax=axes[0,0], zorder=1)
    axes[0,0].set_title("Simplicial complex with original signal")
    plt.colorbar(axes[0,0].collections[1], ax=axes[0,0])
    axes[0,0].set_xticks([])
    axes[0,0].set_yticks([])
    
    
    plot_nodes(s0, points,ax=axes[1,0], zorder=3,s=size_nodes)
    plot_edges(s1phi,points,X, ax=axes[1,0], zorder=2,linewidths=size_lines)
    if len(X)>2: plot_triangles_plain(color_tri,points,X, ax=axes[1,0], zorder=1)
    axes[1,0].set_title("Simplicial complex with reconstructed signal")
    cbar=plt.colorbar(axes[1,0].collections[0], ax=axes[1,0])
    #cbar.set_ticks([]) 
    #cbar.set_ticklabels([])
    axes[1,0].set_xticks([])
    axes[1,0].set_yticks([])
    
    
    
    plot_nodes(s0, points,ax=axes[1,1], zorder=3,s=size_nodes)
    plot_edges_normed(s1loss,points,X, ax=axes[1,1], zorder=2,linewidths=size_lines)
    if len(X)>2: plot_triangles_plain(color_tri,points,X, ax=axes[1,1], zorder=1)
    axes[1,1].set_title("Signal Loss")
    plt.colorbar(axes[1,1].collections[0], ax=axes[1,1])
    axes[1,1].set_xticks([])
    axes[1,1].set_yticks([])
        
    down,up=from_indices_to_color(X,collapses,dimq)
    if type_collapse=="up":
        s1=np.array([mcolors.to_rgba('black')]*len(X[1]))
        s1[down]=mcolors.to_rgba('firebrick')
        if len(X)>2: 
            s2=np.array([mcolors.to_rgba(color_tri)]*len(X[2]))
            s2[up]=mcolors.to_rgba('darksalmon')
    if type_collapse=="down":
        s0=np.array([mcolors.to_rgba('black')]*len(X[0]))
        s0[down]=mcolors.to_rgba('firebrick')
        s1=np.array([mcolors.to_rgba("black")]*len(X[1]))
        s1[up]=mcolors.to_rgba('darksalmon') #sandybrown
    plot_nodes(s0, points,ax=axes[2,0], zorder=3,s=size_nodes)
    plot_edges_plain(s1,points,X, ax=axes[2,0], zorder=2,linewidths=size_lines)
    if len(X)>2: plot_triangles_plain(s2,points,X, ax=axes[2,0], zorder=1)
    axes[2,0].set_title("Collapsed Simplices")
    plt.colorbar(axes[2,0].collections[0], ax=axes[2,0])
    axes[2,0].set_xticks([])
    axes[2,0].set_yticks([])
    
    
    
    cp=list(collapsed_X[0].keys())
    cp=set([list(c)[0] for c in cp])
    tot_p=list(X[0].keys())
    tot_p=set([list(c)[0] for c in tot_p])
    rp=list(tot_p.intersection(cp))
    plot_nodes(np.array(s0)[rp], np.array(points)[rp],ax=axes[0,1], zorder=3,s=size_nodes)
    plot_edges_normed(s1c,points,collapsed_X, ax=axes[0,1], zorder=2,linewidths=size_lines)
    if len(X)>2: plot_triangles_plain(color_tri,points,collapsed_X, ax=axes[0,1], zorder=1)
    axes[0,1].set_title("Collapsed complex")
    plt.colorbar(axes[0,1].collections[0], ax=axes[0,1])
    axes[0,1].set_xticks([])
    axes[0,1].set_yticks([])
        
    axes[-1,-1].set_visible(False)
    
    
    
def plot_optimal_edgecollapses(X,dimq,signal,steps,points,color_tri,size_nodes,size_lines):
    s0 = ['black']*len(X[0])#np.zeros(len(simplices[0]))
    s1 = signal
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20, 15))
    plot_nodes(s0, points,ax=axes[0,0], zorder=3,s=size_nodes)
    plot_edges(s1,points,X, ax=axes[0,0], zorder=2,linewidths=size_lines)
    plot_triangles_plain(color_tri,points,X, ax=axes[0,0], zorder=1)
    axes[0,0].set_title("Simplicial complex with original signal")
    plt.colorbar(axes[0,0].collections[0], ax=axes[0,0])
    axes[0,0].set_xticks([])
    axes[0,0].set_yticks([])
    

    all_X,collapses,all_losses,total_loss,all_signals,s1=sequence_optimal_collapses(X,dimq,signal,steps)
    plot_nodes(s0, points,ax=axes[0,1], zorder=3,s=size_nodes)
    plot_edges(s1,points,X, ax=axes[0,1], zorder=2,linewidths=size_lines)
    plot_triangles_plain(color_tri,points,X, ax=axes[0,1], zorder=1)
    axes[0,1].set_title("Simplicial complex with reconstructed signal")
    plt.colorbar(axes[0,1].collections[0], ax=axes[0,1])
    axes[0,1].set_xticks([])
    axes[0,1].set_yticks([])
        
    signal_loss=np.abs(signal-s1)
    plot_nodes(s0, points,ax=axes[1,0], zorder=3,s=size_nodes)
    plot_edges(signal_loss,points,X, ax=axes[1,0], zorder=2,linewidths=size_lines)
    plot_triangles_plain(color_tri,points,X, ax=axes[1,0], zorder=1)
    axes[1,0].set_title("Signal Loss")
    plt.colorbar(axes[1,0].collections[0], ax=axes[1,0])
    axes[1,0].set_xticks([])
    axes[1,0].set_yticks([])
        
    down,up=from_indices_to_color(X,collapses,dimq)
    s1=np.array([mcolors.to_rgba('black')]*len(X[1]))
    s1[down]=mcolors.to_rgba('firebrick')
    s2=np.array([mcolors.to_rgba(color_tri)]*len(X[2]))
    s2[up]=mcolors.to_rgba('sandybrown')
    plot_nodes(s0, points,ax=axes[1,1], zorder=3,s=size_nodes)
    plot_edges_plain(s1,points,X, ax=axes[1,1], zorder=2,linewidths=size_lines)
    plot_triangles_plain(s2,points,X, ax=axes[1,1], zorder=1)
    axes[1,1].set_title("Collapsed Simplices")
    plt.colorbar(axes[1,1].collections[0], ax=axes[1,1])
    axes[1,1].set_xticks([])
    axes[1,1].set_yticks([])
    
def plot_hodge_decomp(X,s1,kX,phispsis,trange,type_collapse='up',c1='teal',c2='coral'):        
    
    boundaries=kX
    ups,downs,laplacians=build_up_down_laplacians(boundaries)
    
    
    down=downs[1]
    lap=laplacians[1]
    vh,vech=sparse.linalg.eigsh(lap, 1, which='SM')
    basis_h=vech[:,np.where(vh<10**(-6))[0]]


    if len(ups)>1:
        up=ups[1]
        vup,vecup=np.linalg.eigh(up.toarray(), UPLO='L')
        if trange==None:
            uprange=len(np.where(vup>10**(-6))[0])
        if trange!=None:
            uprange=trange
        basis_up=vecup[:,np.where(vup>10**(-6))[0]][:,:uprange]
    
    vdown,vecdown=np.linalg.eigh(down.toarray(), UPLO='L')
    if trange==None:
    
        downrange=len(np.where(vdown>10**(-6))[0])
    if trange!=None:
        
        downrange=trange
        
    
    
    
    basis_down=vecdown[:,np.where(vdown>10**(-6))[0]][:,:downrange]
    

    hodge_basis=basis_h
    if type_collapse=="up":
        hodge_basis=np.hstack((basis_h,basis_down))
        if len(ups)>1:
            hodge_basis=np.hstack((hodge_basis,basis_up))

    if type_collapse=="down":
        if len(ups)>1:
            hodge_basis=np.hstack((basis_h,basis_up))
        hodge_basis=np.hstack((hodge_basis,basis_down))

    h_dec=hodge_basis.T@s1
    h_dec_reconstruction=hodge_basis.T@(phispsis)

    fig = plt.figure(figsize=(11,5))
    l=len(h_dec)
    w=l/(l*2)
    x=np.arange(l)
    plt.bar(x, h_dec, width=w,label="Original signal",color=c1)
    x1=x+w
    plt.bar(x1, h_dec_reconstruction, width=w,label="Reconstructed signal",color=c2)
    plt.xticks([], [])
    plt.title("Projection of the signal on the Hodge decomposition")
    plt.legend()
    
    
def height_function(X,points):
    lines=get_positions(points,X,1)
    height_fun=[]
    for i in range(len(lines)):
        height_fun.append((lines[i][0][1]+lines[i][1][1])/2)
    return height_fun



def compute_min_max(signal,phipsis,collapsed_signal):
    s1 = signal
    s1phi=phipsis
    s1loss=np.abs(s1-s1phi)
    s1c=collapsed_signal
    a=[s1,s1phi,s1loss,s1c]
    tot_max=max(map(max, a))
    tot_min=min(map(min, a))
    return tot_max,tot_min




   