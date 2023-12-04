import random
import matplotlib.pyplot as plt

def k_clusters_with_l_points(k:int, l:int, rangeX, rangeY, variance):
    points = list()
    for i in range(k):
        x = random.random()*rangeX
        y = random.random()*rangeY
        points.append((x,y))
        for j in range(l):
            a = x + (random.random()-0.5)*variance
            b = y + (random.random()-0.5)*variance
            points.append((a,b))
    return points

rangeX = 500
rangeY = 500
variance = 40

output = k_clusters_with_l_points(5,20,rangeX,rangeY,variance)

print(output) 


"""x_axis = list()
y_axis = list()
for i in output:
    x_axis.append(i[0])
    y_axis.append(i[1])
plt.scatter(x_axis,y_axis)
plt.show()"""