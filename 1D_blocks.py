import cmath
import math
import matplotlib.pyplot as plt
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
    def set_antennas(self, I):
        self.antennas[0].set_values(I[0])
        self.antennas[1].set_values(I[1])
    def get_antenna_vals(self):
        return [self.antennas[i].get_phase_amp() for i in range(2)]

class BlockArray:
    def __init__(self, ps, t):
        self.target = t
        self.blocks = [Block(p) for p in ps]
        self.num_blocks = len(self.blocks)
        a_pos = [0] * (2 * self.num_blocks)
        for b in range(self.num_blocks):
            center = self.blocks[b].get_pos()
            offset = self.blocks[b].get_spacing() / 2
            a_pos[b * 2] = center - offset
            a_pos[(b * 2) + 1] = center + offset
        self.a_pos = a_pos

    def get_array_factor(self, theta):
        sum = 0
        I_0 = self.blocks[0].get_antennas()[0].get_excitation()
        for a in range(self.num_blocks * 2):
            I_n = self.blocks[int(a / 2)].get_antennas()[a % 2].get_excitation()
            sum += (I_n / I_0) * cmath.exp(1j * k * self.a_pos[a] * math.cos(theta))
        return sum

    def get_solution(self):
        for b in range(self.num_blocks):
            I = [cmath.exp(1j * -k * self.a_pos[(2 * b) + i]
                * math.cos(self.target)) for i in range(2)]
            self.blocks[b].set_antennas(I)

    def output_results(self):
        # amps = [0] * self.num_blocks * 2
        # phases = [0] * self.num_blocks * 2
        # for b in range(0, self.num_blocks * 2, 2):
        #     a1, a2 = self.blocks[int(b / 2)].get_antenna_vals()
        #     phases[b], amps[b] = a1
        #     phases[b + 1], amps[b + 1] = a2
        #
        # block_nums = list(sum(zip([str(i + 1) for i in range(self.num_blocks)], [""] * self.num_blocks), ())[:-1]) + ['']
        # # output antenna configurations
        # print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
        # array_config = zip(*[block_nums, amps, phases])
        # print (tabulate(array_config, headers=['Block #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

        amps = [0] * self.num_blocks
        phases = [0] * self.num_blocks
        print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
        for b in range(self.num_blocks):
            a1, a2 = self.blocks[b].get_antenna_vals()
            phase1, amp1 = a1
            phase2, amp2 = a2
            amps[b] = (amp1, amp2)
            phases[b] = (phase1, phase2)
            print ('Block ' + str(b + 1) + ':')
            array_config = zip(*[[amp1, amp2], [phase1, phase2]])
            print (tabulate(array_config, headers=['Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

        # # get array factor at each target angle
        # target_results = []
        # for target in targets:
        #     target_results.append(abs(get_array_factor(math.radians(target))))
        # closest_peaks = []
        # for angle in closest:
        #     closest_peaks.append(abs(get_array_factor(math.radians(angle))))
        #
        # # output results
        # print ('\n\033[1m Resulting Beamform:\033[0m\n')
        # print ('--- Best uniform peak array factor: ' + str(M / T))
        # print ('--- Runtime: ' + str(runtime) + ' seconds\n')
        # results = zip(*[[x + 1 for x in range(T)], targets, closest, target_results, closest_peaks])
        # print (tabulate(results, headers=['Target #', 'Target Angle', 'Closest Sample', 'Array Factor', 'AF at Closest'], tablefmt='orgtbl') + '\n')
        # diffs = [abs(targets[i] - closest[i]) for i in range(T)]
        # print ('--- Mean difference target vs. closest: ' + str(sum(diffs)/T))
        # print ('--- Max difference target vs. closest: ' + str(max(diffs)))
        # print ('--- Mean array factor at target: ' + str(sum(target_results)/T))
        # print ('--- Max array factor at target: ' + str(max(target_results)))
        # print ('--- Min array factor at target: ' + str(min(target_results)) + '\n')

    def visualize(self, show_block_positions=False):
        if show_block_positions:
            plt.scatter(self.a_pos, [0] * (self.num_blocks * 2))
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



# array1 = BlockArray([-20, -lmbda, -100, -90, -75, -60, 1.3 * lmbda, 30, 35, 40, 50], math.radians(125))
array1 = BlockArray([-0, -lmbda, lmbda, -2 * lmbda, 2 * lmbda, -3 * lmbda, 3 * lmbda, -4*lmbda, 4*lmbda, -5*lmbda, 5*lmbda], math.radians(125))
array1.get_solution()
array1.output_results()
array1.visualize(True)
