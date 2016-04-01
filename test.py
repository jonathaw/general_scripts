from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# X, Y, Z = axes3d.get_test_data(0.05)
#
# print('X', X, type(X))
# print('Y', Y)
# print('Z', Z)


X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
Y = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
Z = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
# plt.zlabel('Z')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
