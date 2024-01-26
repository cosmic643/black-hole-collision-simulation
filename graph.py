import matplotlib.pyplot as plt

x = []
y = []
with open('time.txt','r') as x_axis:
    for i in range(1000):
        x.append(float(x_axis.readline()));

with open('pno.txt','r') as y_axis:
    for i in range(1000):
        y.append(float(y_axis.readline()));
    
print(x)
print(y)

plt.plot(x,y,'g')
plt.show()

