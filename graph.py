import matplotlib.pyplot as plt

x = []
y = []

with open('time.txt','r') as temp:
    l = temp.read()
    l_t = l.split('\n')

data_size = len(l_t)
print(data_size)

with open('time.txt','r') as x_axis:
    for i in range(data_size-1):
        x.append(float(x_axis.readline()))

with open('pno.txt','r') as y_axis:
    for i in range(data_size-1):
        y.append(float(y_axis.readline()))

with open('black_hole_system_info.txt','r') as system_data:
    mean_mass = float(system_data.readline())
    sd_mass = float(system_data.readline())
    mean_pos = float(system_data.readline())
    sd_pos = float(system_data.readline())
    n = float(system_data.readline())
    itr = float(system_data.readline())
    TIME_SCALE = float(system_data.readline())
    merger_radius = float(system_data.readline())

mean_mass = 'Mean Mass: ' + str(mean_mass)
sd_mass = 'SD of Mass: ' + str(sd_mass)
mean_pos = 'Mean position: ' + str(mean_pos)
sd_pos = 'SD of position: ' + str(sd_pos)
n = str(int(n))
itr = str(int(itr))
TIME_SCALE = "Time per Iteration: " + str(TIME_SCALE)
merger_radius = 'Merging Radius: ' + str(merger_radius)

data = mean_mass + '\n' + sd_mass + '\n' + mean_pos + '\n' + sd_pos + '\n' + TIME_SCALE + '\n' + merger_radius

plt.xlabel("Iterations (Total iterations: " + itr + ')', fontsize = 10)
plt.ylabel("Number of Particles (Initially: " + n + ')',fontsize = 10)
plt.text(400,60, data, fontsize = 12,bbox = dict(facecolor = 'white', alpha = 0.5))

plt.plot(x,y,'g')
plt.show()

