#!/usr/bin/env python3.5
import numpy as np
import matplotlib.pyplot as plt
import random

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


x = [random.randint(1, 4) for a in range(500)]
print(x)
col_dict = {1: 'black', 2: 'red', 3: 'green', 4: 'purple'}
for i, x_ in enumerate(x):
    plt.hlines(1, i, i+1, colors=col_dict[x_], lw=20)
plt.plot(range(500), x)

plt.show()