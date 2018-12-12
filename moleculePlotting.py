'''
code helps visualize the differences between the original molecule (in blue) against the new molecule (in red)
'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from errorMinSubRoutine import solution
from targetOutputBase import originalMatrix

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

x1data = solution[1,:]
y1data = solution[2,:]
z1data = solution[3,:]

x2data, y2data, z2data = originalMatrix[1,:],originalMatrix[2,:],originalMatrix[3,:]

ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
ax.set_zlabel('z axis')

ax.scatter(x1data, y1data, z1data, c = 'r')
ax.scatter(x2data, y2data, z2data, c = 'b')

plt.show()
