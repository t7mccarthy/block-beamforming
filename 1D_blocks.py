import cmath
import math
import matplotlib.pyplot as plt

lmbda = 5.168835482759
wave_num = 2 * math.pi / lmbda

class Antenna:
    def __init__(self):
        self.phase = None
        self.amplitude = None
        self.excitation = 1
    def set_values(self, I):
        self.excitation = I
        x = I.real
        y = I.imag
        self.phase = math.atan2(y, x)
        self.amplitude = x / math.cos(self.phase)
    def get_phase_amp(self):
        return (self.phase, self.amplitude)
    def get_excitation(self):
        return self.excitation

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
    def set_antennas(self, a1, a2):
        self.antennas[0].set_values(a1)
        self.antennas[1].set_values(a2)

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

    def get_array_factor(self, theta):
        sum = 0
        n = -(self.num_blocks * 2 - 1)/2
        I_0 = self.blocks[0].get_antennas()[0].get_excitation()
        for b in range(self.num_blocks):
            for i in range(2):
                a = 2 * b + i
                d = self.antenna_pos[a]/n
                I_n = self.blocks[b].get_antennas()[i].get_excitation()
                sum += (I_n / I_0) * cmath.exp(1j * n * wave_num * d * math.cos(theta))
                n += 1
        return sum

    def get_solution(self):
        n = -(self.num_blocks * 2 - 1)/2
        for b in range(self.num_blocks):
            I = [0] * 2
            for i in range(2):
                n_i = n + i
                d = self.antenna_pos[(2 * b) + i]/n_i
                a = wave_num * d * math.cos(self.target)
                I[i] = cmath.exp(1j * -n_i * a)
            n += 2
            self.blocks[b].set_antennas(I[0], I[1])

    def visualize(self, show_block_positions=False):
        if show_block_positions:
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



# array1 = BlockArray([-20, -lmbda, 1.3 * lmbda, 30, 35, 40, 50], math.radians(125))
array1 = BlockArray([-20, -lmbda, 1.3 * lmbda, 30, 35, 40, 50], math.radians(125))
array1.get_solution()
array1.visualize()
