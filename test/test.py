from abc import ABC, abstractmethod  # abstract base classes (ABC)

import matplotlib.pyplot as plt
import numpy as np


class Point():
    def __init__(self, xy=None):
        if xy is None:
            self.x, self.y = np.random.random(2)
        else:
            self.x, self.y = xy

    def __repr__(self):
        return "Point: %.2f, %.2f" % (self.x, self.y)

    def distance(self, point):
        dx = (self.x - point.x)
        dy = (self.y - point.y)
        dist = (dx**2 + dy**2) ** .5
        return dist


def estimate_pi(n_samples):
    pt_center = Point((.5, .5))
    pts = [Point() for _ in range(n_samples)]
    dists = [pt.distance(pt_center) for pt in pts]
    ratio = np.mean(np.array(dists) <= .5)
    pi = ratio / (.5**2)
    return pi


outs = np.zeros((100, 15))
for i in range(15):
    outs[:, i] = [estimate_pi(2**i) for _ in range(100)]

plt.boxplot(outs)
plt.axhline(y=3.141596)


class Creatures(ABC):
    @abstractmethod
    def legs(self):
        pass

    @abstractmethod
    def can_fly(self):
        pass


class Bird(Creatures):
    def legs(self):
        return 2

    def can_fly(self):
        return True

    def jump(self):
        return 3


class Human(Creatures):

    def __init__(self, name):
        self.

    def has_legs(self):
        return 2

    def can_fly(self):
        return True

    def jump(self):
        return 10

    def jump(self, use_tools=True):
        if use_tools:

        return 10


b = Bird()


b.legs()


class Polygon(ABC):

    @abstractmethod
    def noofsides(self):
        pass


class Triangle(Polygon):

    # overriding abstract method
    def noofsides(self):
        print("I have 3 sides")
