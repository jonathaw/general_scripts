#!/usr/bin/env python3.5
"""
Animation of Elastic collisions with Gravity

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

class Particle:
    def __init__(self, loc, speed, col):
        self.loc = loc
        self.speed = speed
        self.col = col

    def step(self):
        self.loc[0] += self.speed[0]
        self.loc[1] += self.speed[1]


class Box:
    def __init__(self, n):
        self.particles = []
        for i in range(n):
            self.particles.append(Particle([0, 0], [1, 1], 'k'))

    def step(self):
        for particle in self.particles:
            particle.step()

def init():
    global particles
    particles.set_data([], [])
    return particles,


# animation function.  This is called sequentially
def animate(i):
    global particles, box
    box.step()
    return particles,

fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))

box = Box(10)
particles, = ax.plot([], [])

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)
plt.show()
