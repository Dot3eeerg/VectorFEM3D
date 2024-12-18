import matplotlib.pyplot as plt
import numpy as np

x, y, Az, u, b, s, g = [], [], [], [], [], [], []

#with open("time") as file:
#    for line in file:
#        x.append(float(line))

with open("timeNew") as file:
    for line in file:
        x.append(float(line))

with open("x") as file:
    for line in file:
        u.append(float(line))

with open("y") as file:
    for line in file:
        b.append(float(line))

with open("r0") as file:
    for line in file:
        u.append(float(line))


x = list(dict.fromkeys(x))

# plot
fig, ax = plt.subplots()

#X, Y = np.meshgrid(x, u)
ax.plot(x, u, marker=5, label='A phi')
ax.plot(x, s, label='A vect')

ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()
