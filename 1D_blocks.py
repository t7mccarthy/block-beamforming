import cmath
import math
import matplotlib.pyplot as plt
from statistics import mean
from tabulate import tabulate


lmbda = 5.168835482759
k = 2 * math.pi / lmbda

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
        return (math.degrees(self.phase), self.amplitude)

class Block:
    def __init__(self, p, s = lmbda / 2):
        self.pos = p
        self.spacing = s
        self.antennas = (Antenna(), Antenna())
    def set_antennas(self, I):
        self.antennas[0].set_values(I[0])
        self.antennas[1].set_values(I[1])
    def get_antenna_vals(self):
        return ([self.antennas[0].amplitude, self.antennas[1].amplitude], [self.antennas[0].phase, self.antennas[1].phase])
    def get_excitation(self, i):
        return self.antennas[i].excitation

class BlockArray:
    def __init__(self, ps, t):
        self.target = t
        self.blocks = [Block(p) for p in ps]
        self.num_blocks = len(self.blocks)
        a_pos = [0] * (2 * self.num_blocks)
        for b in range(self.num_blocks):
            center = self.blocks[b].pos
            offset = self.blocks[b].spacing / 2
            a_pos[b * 2] = center - offset
            a_pos[(b * 2) + 1] = center + offset
        self.a_pos = a_pos
        self.afs = False

    def get_array_factor(self, theta):
        sum = 0
        I_0 = self.blocks[0].get_excitation(0)
        for a in range(self.num_blocks * 2):
            I_n = self.blocks[int(a / 2)].get_excitation(a % 2)
            sum += (I_n / I_0) * cmath.exp(1j * k * self.a_pos[a] * math.cos(theta))
        return sum

    def calc_afs(self):
        if not self.afs:
            y = [0] * 181
            for i in range(181):
                y[i] = abs(self.get_array_factor(math.radians(i)))
            self.afs = y[:]

    def solve(self):
        for b in range(self.num_blocks):
            I = [cmath.exp(1j * -k * self.a_pos[(2 * b) + i] * math.cos(self.target)) for i in range(2)]
            self.blocks[b].set_antennas(I)

    def output_results(self):
        print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
        for b in range(self.num_blocks):
            amps, phases = self.blocks[b].get_antenna_vals()
            print ('Block ' + str(b + 1) + ':')
            array_config = zip(*[[1, 2], amps, phases])
            print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')
        print ('\n\033[1m Resulting Beamform:\033[0m\n')
        print ('--- Best expected array factor: ' + str(self.num_blocks * 2))
        print ('--- Array factor at target: ' + str(abs(self.get_array_factor(self.target))))
        self.calc_afs()
        print ('--- Mean array factor: ' + str(mean(self.afs)))
        print ('\nFINISHED\n')

    def visualize(self, show_block_positions=False):
        if show_block_positions:
            plt.scatter(self.a_pos, [0] * (self.num_blocks * 2))
            plt.title('Antenna Positions')
            plt.xlabel('Position')
            plt.show()
        self.calc_afs()
        plt.plot(range(181), self.afs)
        plt.title('Block-Generated Beamform')
        plt.xlabel('Angle from Array (Degrees)')
        plt.ylabel('Array Factor')
        plt.show()



array1 = BlockArray([-20, -lmbda, -100, -90, -75, -60, 1.3 * lmbda, 30, 35, 40, 50], math.radians(125))
# array1 = BlockArray([-0, -lmbda, lmbda, -2 * lmbda, 2 * lmbda, -3 * lmbda, 3 * lmbda, -4*lmbda, 4*lmbda, -5*lmbda, 5*lmbda], math.radians(125))
array1.solve()
array1.output_results()
array1.visualize(True)
