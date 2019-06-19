import matplotlib.pyplot as plt

lmbda = 5.168835482759

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

    def visualize(self):
        plt.scatter(self.antenna_pos, [0] * (self.num_blocks * 2))
        plt.title('Antenna Positions')
        plt.xlabel('Position')
        plt.show()

array1 = BlockArray([-10, 20, 28, 40, 70], 10)
array1.visualize()
