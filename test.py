#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import interpolate

def MakeHydrophobicityGrade():
    """
    :return: returns a dictionary of the polynom values for each residue
    """
    global hydrophobicity_polyval
    try:
        hydrophobicity_grade = open(ELAZAR_POLYVAL_PATH, 'r')
    except:
        hydrophobicity_grade = open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/mother_fucker_MET_LUE_VAL_sym_k_neg.txt', 'r')
    hydrophobicity_polyval = {}
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('click on points')
#
# line, = ax.plot(np.random.rand(100), 'o', picker=5)  # 5 points tolerance
#
# def onpick(event):
#     thisline = event.artist
#     xdata = thisline.get_xdata()
#     ydata = thisline.get_ydata()
#     ind = event.ind
#     points = tuple(zip(xdata[ind], ydata[ind]))
#     print('onpick points:', points)
#
# fig.canvas.mpl_connect('pick_event', onpick)
#
# plt.show()


# x = [random.randint(1, 4) for a in range(500)]
# print(x)
# col_dict = {1: 'black', 2: 'red', 3: 'green', 4: 'purple'}
# for i, x_ in enumerate(x):
#     plt.hlines(1, i, i+1, colors=col_dict[x_], lw=20)
# plt.plot(range(500), x)
#
# plt.show()

# x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
# y = np.sin(x)
# tck = interpolate.splrep(x, y, s=0)
# xnew = np.arange(0, 2*np.pi, np.pi/50)
# ynew = interpolate.splev(xnew, tck, der=0)
profiles = MakeHydrophobicityGrade()
x = np.arange(-50, +50, 0.1)
Z = np.arange(-15, 15, 0.1)


# print(x, len(x))
# print(y, len(y))

plt.figure()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
for i, aa in enumerate(list('ACDEFGHIKLMNPQRSTVWY')):
    y = [0] * (35 * 10) + [np.polyval(profiles[aa], z) for z in Z] + [0] * (35 * 10)
    if aa not in ['V', 'M', 'H']:
        tck = interpolate.splrep(x, y, s=25)
    else:
        tck = interpolate.splrep(x, y, s=5)
    print('for %s tck is:' % (aa))
    print('t', len(tck[0]))
    print('c', len(tck[1]))
    print('k', tck[2])
    xnew = np.arange(-50, 50, 0.05)
    ynew = interpolate.splev(xnew, tck, der=0)
    plt.subplot(5, 4, 1+i)
    plt.plot(x, y, c='k')
    plt.plot(xnew, ynew, c='r')
    plt.vlines(-15, -2, 3, color='grey', linestyles='dashed')
    plt.vlines(15, -2, 3, color='grey', linestyles='dashed')
    plt.title(aa.upper())
    plt.ylim([-2, 3])


# plt.figure()
# plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
# plt.legend(['Linear', 'Cubic Spline', 'True'])
# plt.axis([-0.05, 6.33, -1.05, 1.05])
# plt.title('Cubic-spline interpolation')
plt.show()
