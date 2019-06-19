import cmath
import math
import matplotlib.pyplot as plt

lmbda = 5.168835482759
wave_num = 2 * math.pi / lmbda

class Antenna:
    def __init__(self):
        self.phase = None
        self.amplitude = None
    def set_values(self, p, a):
        self.phase = p
        self.amplitude = a
    def get_values(self):
        return (self.phase, self.amplitude)

class Block:
    def __init__(self, p, s = lmbda / 2):
        self.pos = p
        self.spacing = s
        self.antennas = (Antenna(), Antenna())
    def get_pos(self):
        return self.pos
    def get_spacing(self):
        return self.spacing
    def get_antennas(self):
        return self.antennas

class BlockArray:
    def __init__(self, ps, t):
        self.target = t
        self.blocks = [Block(p) for p in ps]
        self.num_blocks = len(self.blocks)
        antenna_pos = [0] * (2 * self.num_blocks)
        for b in range(self.num_blocks):
            center = self.blocks[b].get_pos()
            offset = self.blocks[b].get_spacing() / 2
            antenna_pos[(b * 2) - 1] = center - offset
            antenna_pos[b * 2] = center + offset
        self.antenna_pos = antenna_pos
        self.phases = [0] * (self.num_blocks * 2)

    def change_target(self, t):
        self.target = t

    def get_array_factor(self, theta):
        # d = lmbda / 2
        sum = 0
        n = -(self.num_blocks * 2 - 1)/2
        for a in range(self.num_blocks * 2):
            # d = lmbda/2
            d = self.antenna_pos[a]/n
            r = wave_num * d * math.cos(self.target)
            sum += cmath.exp(1j * n * (wave_num * d * math.cos(theta) - r))
            n += 1
        return sum

    # def get_solution(self):
    #     d = lmbda / 2
    #     a = k * d * math.cos(self.target)
    #

    # def find_best(self):


    def visualize(self):
        plt.scatter(self.antenna_pos, [0] * (self.num_blocks * 2))
        plt.title('Antenna Positions')
        plt.xlabel('Position')
        plt.show()

        y = [0] * 181
        for i in range(181):
            y[i] = abs(self.get_array_factor(math.radians(i)))
        plt.plot(range(181), y)
        plt.title('Block-Generated Beamform')
        plt.xlabel('Angle from Array (Degrees)')
        plt.ylabel('Array Factor')
        plt.show()



array1 = BlockArray([-20, -lmbda, 1.3 * lmbda, 30, 35, 40], math.radians(25))
# array1.find_best()
array1.visualize()
