import pyhepmc

import matplotlib.pyplot as plt

f = pyhepmc.io.open("DH_100k.hepmc", mode="r")

x = []
y = []
z = []

for i_event, event in enumerate(f):
    if not event:
        continue
    print(i_event)
    if i_event > 631:
        break
    if event.vertices:
        x.append(event.vertices[0].position.x)
        y.append(event.vertices[0].position.y)
        z.append(event.vertices[0].position.z)
        print(i_event, " done")


plt.figure()
plt.hist(z, bins = 100)
plt.xlabel("z")

plt.figure()
plt.hist(x, bins = 100)
plt.xlabel("x")

plt.figure()
plt.hist(y, bins = 100)
plt.xlabel("y")

plt.figure()
plt.hist2d(x, y, bins = 100)
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar()

plt.show()
