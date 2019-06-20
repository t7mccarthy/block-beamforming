import cmath
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
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
    def __init__(self, p, x_s = lmbda / 2, y_s = lmbda / 2):
        self.pos = p
        self.x_spacing = x_s
        self.y_spacing = y_s
        self.antennas = (Antenna(), Antenna(), Antenna(), Antenna())
    def set_antennas(self, I):
        for i in range(4):
            self.antennas[i].set_values(I[i])
    def get_antenna_vals(self):
        return ([(a.phase, a.amplitude) for a in self.antennas])
    def get_excitation(self, i):
        return self.antennas[i].excitation

class BlockArray:
    def __init__(self, ps, t, p):
        self.target = (t, p)
        self.blocks = [Block(p) for p in ps]
        self.num_blocks = len(self.blocks)
        a_pos = [0] * (4 * self.num_blocks)
        for b in range(self.num_blocks):
            x, y = self.blocks[b].pos
            x_offset = self.blocks[b].x_spacing / 2
            y_offset = self.blocks[b].y_spacing / 2
            a_pos[b * 4] = (x - x_offset, y + y_offset)
            a_pos[(b * 4) + 1] = (x + x_offset, y + y_offset)
            a_pos[(b * 4) + 2] = (x - x_offset, y - y_offset)
            a_pos[(b * 4) + 3] = (x + x_offset, y - y_offset)
        self.a_pos = a_pos
        self.afs = False

    def get_array_factor(self, theta, phi):
        sum = 0
        I_0 = self.blocks[0].get_excitation(0)
        for a in range(self.num_blocks * 2):
            I_n = self.blocks[int(a / 2)].get_excitation(a % 2)
            sum += (I_n / I_0) * cmath.exp(1j * k * self.a_pos[a] * math.cos(theta))
        return sum


    def get_array_factor (self, theta, phi):
        sum = 0
        I_00 = self.blocks[0].get_excitation(0)
        for a in range(self.num_blocks * 4):
            I_uv = self.blocks[int(a / 4)].get_excitation(a % 4)
            sum += (I_uv/I_00) * cmath.exp(1j * k * math.sin(theta) * (self.a_pos[a][0] * math.cos(phi) + self.a_pos[a][1] * math.sin(phi)))
        return abs(sum)

    def calc_afs(self):
        if not self.afs:
            afs = [[0] * 360 for i in range(90)]
            for t in range(90):
                for p in range(360):
                    y[i] = abs(self.get_array_factor(math.radians(t), math.radians(p)))
            self.afs = afs[:]

    # def solve(self):
    #     for b in range(self.num_blocks):
    #         I = [cmath.exp(1j * -k * self.a_pos[(2 * b) + i] * math.cos(self.target)) for i in range(2)]
    #         self.blocks[b].set_antennas(I)

    # def output_results(self):
    #     print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
    #     for b in range(self.num_blocks):
    #         amps, phases = self.blocks[b].get_antenna_vals()
    #         print ('Block ' + str(b + 1) + ':')
    #         array_config = zip(*[[1, 2], amps, phases])
    #         print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')
    #     print ('\n\033[1m Resulting Beamform:\033[0m\n')
    #     print ('--- Best expected array factor: ' + str(self.num_blocks * 2))
    #     print ('--- Array factor at target: ' + str(abs(self.get_array_factor(self.target))))
    #     self.calc_afs()
    #     print ('--- Mean array factor: ' + str(mean(self.afs)))
    #     print ('\nFINISHED\n')

    def visualize(self, show_block_positions=False):
        if show_block_positions:
            plt.scatter([x[0] for x in self.a_pos], [x[1] for x in self.a_pos])
            plt.title('Antenna Positions')
            plt.xlabel('X Position')
            plt.ylabel('Y Position')
            plt.show()

        dim = 100
        theta, phi = np.linspace(0, np.pi, dim), np.linspace(0, 2 * np.pi, dim)
        THETA, PHI = np.meshgrid(theta, phi)
        r = [[0]*dim for _ in range(dim)]
        for i in range(dim):
            for j in range(dim):
                r[i][j]=self.get_array_factor(THETA[i][j], PHI[i][j])
        R = np.array(r)
        X = R * np.sin(THETA) * np.cos(PHI)
        Y = R * np.sin(THETA) * np.sin(PHI)
        Z = -abs(R * np.cos(THETA))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        plot = ax.plot_surface(X, Y, Z)

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.show()
        # self.calc_afs()
        # plt.plot(range(181), self.afs)
        # plt.title('Block-Generated Beamform')
        # plt.xlabel('Angle from Array (Degrees)')
        # plt.ylabel('Array Factor')
        # plt.show()



# array1 = BlockArray([(0, 0), (10, 8), (-17, 20), (-7, 0), (4, -25)], math.radians(30), math.radians(120))
array1 = BlockArray([(lmbda/2, lmbda/2), (-lmbda/2, -lmbda/2), (lmbda/2, -lmbda/2), (-lmbda/2, lmbda/2)], math.radians(30), math.radians(120))
# array1.solve()
# array1.output_results()
array1.visualize(True)
