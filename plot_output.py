#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import os
import csv

f_list = []
for f in os.listdir('./'):#find all data files that starts with, say "forces_b2".
    if f.startswith('output_'):
        f_list.append(f)
        f_list.sort()
print "Found these output files:\n"
print f_list
print '\n'

for f in f_list:
    print "reading: ",f
    with open(f, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        history = list(reader)
        history = np.array(history)
        history = history.astype(np.float)#convert string to floats

    fig, axes = plt.subplots(2,1, sharex=True, figsize=(12, 6))
    fig.suptitle('Result of 1D shallow water equation', fontsize=12)
    axes[0].plot(history[:,0], history[:,1],'bo')
    axes[0].legend(['$q^1$'])
    axes[0].set_xlabel("$x$")
    axes[0].set_ylabel("$h$")
    #axes[0].set_title("I like $\pi$")
    margin = 0.1
    y0range = (history[:,1].max()-history[:,1].min()) 
    axes[0].set_ylim(history[:,1].min()-margin*y0range,history[:,1].max()+margin*y0range)
    axes[1].plot(history[:,0], history[:,2],'bo')
    axes[1].legend(['$q^2$'])
    axes[1].set_xlabel("$x$")
    axes[1].set_ylabel("$hu$")
    y1range = (history[:,2].max()-history[:,2].min()) 
    axes[1].set_ylim(history[:,2].min()-margin*y1range,history[:,2].max()+margin*y1range)
    if not('myplot' in os.listdir('./')):
        os.mkdir('myplot')
    plt.savefig('./myplot/'+f+'.png', bbox_inches='tight',dpi=100)
    plt.close()

