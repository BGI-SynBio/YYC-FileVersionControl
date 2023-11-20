import copy
import random
import sys

import math
import numpy
from yyc.utils import log, validity
from yyc.utils.monitor import Monitor
import numpy as np


base_index = {"A": 0, "C": 1, "G": 2, "T": 3}
index_base = {0: "A", 1: "C", 2: "G", 3: "T"}


# noinspection PyProtectedMember
class YYC:
    def __init__(
        self,
        base_reference=None,
        current_code_matrix=None,
        support_bases=None,
        support_spacing=0,
        max_ratio=0.8,
        search_count=100,
        max_homopolymer=math.inf,
        max_content=1,
        min_free_energy=None,
        seed=30
    ):
        """
        The initialization method of YYC.

        :param base_reference: correspondence between base and binary data (RULE 1).
                               Make sure that Two of the bases are 1 and the other two are 0, so there are only 6 case.
        :type: list

        :param current_code_matrix: conversion rule between base and btis based on support and current base (RULE 2).
                                    Label row is the support base, label col is the current base.
                                         A   T   C   G
                                     A   X1  Y1  X2  Y2
                                     T   X3  Y3  X4  Y4
                                     C   X5  Y5  X6  Y6
                                     G   X7  Y7  X8  Y8
                                    Make sure that Xn + Yn = 1 and Xn * Yn = 0, n is in [1, 8].
        :type: list(list)

        :param support_bases: base replenishment before official data.
                              Make sure that the count of support base must more than support spacing.
                              Make sure that the number range of each position is {0, 1, 2, 3}, reference base index.
        :type: string

        :param support_spacing: spacing between support base and current base.
                                When the support base is the front of the current base, the spacing is 0.
        :type: int

        :param max_ratio: the max ratio of 0 or 1.
                          When the (count/length) >= this parameter, we decide that this binary sequence is not good.
        :type: float

        :param max_homopolymer: maximum length of homopolymer.
        :type: int

        :param max_content: maximum content of C and G, which means GC content is in [1 - max_content, max_content].
        :type: float or int

        :param min_free_energy: the free energy of DNA sequence is lower than required min free energy.
        :type: float or none.

        :param seed: random seed.
        :type: int
        """

        # Set default values for Rules 1 and 2 (RULE 495)
        if not base_reference:
            base_reference = [0, 1, 0, 1]
        if not current_code_matrix:
            current_code_matrix = [
                [1, 1, 0, 0],
                [1, 0, 0, 1],
                [1, 1, 0, 0],
                [1, 1, 0, 0],
            ]
        if not support_bases:
            support_bases = [index_base[0]]

        # Assign input data to class variables
        self.base_reference = base_reference
        self.current_code_matrix = current_code_matrix
        self.support_bases = support_bases
        self.support_spacing = support_spacing
        self.max_ratio = max_ratio
        self.search_count = search_count

        self.max_homopolymer = max_homopolymer
        self.max_content = max_content
        self.min_free_energy = min_free_energy

        # Detect parameters correctness
        self._init_check()

        self.file_size = 0
        self.monitor = Monitor()

        self.seed = seed

    def _init_check(self):
        """
        The verification of initialization parameters.

        """
        # Check support bases
        for index in range(len(self.support_bases)):
            if (self.support_bases[index] != "A" and self.support_bases[index] != "T"
                    and self.support_bases[index] != "C" and self.support_bases[index] != "G"):
                log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                           "Only A, T, C, and G can be included as support bases, "
                           "but the support bases[" + str(index) + "] has been detected as "
                           + str(self.support_bases[index] + "!"))

        if len(self.support_bases) < self.support_spacing + 1:
            log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                       "The count of support base needs to be more than support spacing!")

        # Check base reference (rule 1)
        for index in range(len(self.base_reference)):
            if self.base_reference[index] != 0 and self.base_reference[index] != 1:
                log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                           "Only 0 and 1 can be included for base reference, and base_reference[" + str(index)
                           + "] has been detected as " + str(self.base_reference[index]) + "!")
        if sum(self.base_reference) != 2:
            log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                       "Wrong correspondence between base and binary data!")

        positions = []
        for i in range(len(self.base_reference)):
            if self.base_reference[i] == 1:
                positions.append(i)
        for i in range(len(self.base_reference)):
            if self.base_reference[i] != 1:
                positions.append(i)

        # Check current code matrix (rule 2)
        for row in range(len(self.current_code_matrix)):
            for col in range(len(self.current_code_matrix[row])):
                if self.current_code_matrix[row][col] != 0 and self.current_code_matrix[row][col] != 1:
                    log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                               "Only 0 and 1 can be included in the current code matrix, and the current code matrix ["
                               + str(row) + ", " + str(col) + "] has been detected as "
                               + str(self.current_code_matrix[row][col]) + "!")

        for row in range(len(self.current_code_matrix)):
            left = self.current_code_matrix[row][positions[0]]
            right = self.current_code_matrix[row][positions[1]]
            if left + right == 1 and left * right == 0:
                continue
            else:
                log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                           "Wrong current code matrix, the error locations are [" + str(row) + ", " + str(positions[0])
                           + "] and [" + str(row) + ", " + str(positions[1]) + "]! "
                           + "It is required by rule that these two values will have sum of 1 and product of 0.")

            left = self.current_code_matrix[row][positions[2]]
            right = self.current_code_matrix[row][positions[3]]
            if left + right == 1 and left * right == 0:
                continue
            else:
                log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                           "Wrong current code matrix, the error locations are [" + str(row) + ", " + str(positions[2])
                           + "] and [" + str(row) + ", " + str(positions[3]) + "]! "
                           + "It is required by rule that these two values will have sum of 1 and product of 0.")
        # Check max ratio
        if self.max_ratio <= 0.5 or self.max_ratio >= 1:
            log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                       "Wrong max ratio (" + str(self.max_ratio) + ")!")

    # ================================================= encode part ====================================================

    def encode(self, binary_lst, need_log=False):
        """
        Encode DNA sequences from the binary sequences.

        :param binary_lst: generated binary sequences.
                           The element of this list contains only 0 or 1 (non-char).
        :type: list(char)

        :param need_log: show the log.
        :type: bool
        """
        self.monitor.restore()

        if need_log:
            log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                       "Random incorporation and validity testing.")

        data_set = self._adjacent_pairing(binary_lst, need_log)

        self.monitor.restore()

        if need_log:
            log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                       "Convert to DNA sequence string set.")

        dna_sequences = self._synthesis_sequences(data_set, need_log)

        self.monitor.restore()

        return dna_sequences

    def _adjacent_pairing(self, binary_lst, need_log):
        """
        Pair with adjacent binary sequence.

        :param binary_lst: generated binary sequences.
                           The element of this list contains only 0 or 1 (non-char).
        :type: list(char)

        :param need_log: show the log.
        :type: bool
        """
        random.seed(self.seed)

        data_set = []

        total_count = len(binary_lst)

        # index_bit_length = int(len(str(bin(total_count))) - 2)
        index_bit_length = 20

        search_counts = [0 for _ in range(self.search_count + 1)]
        additional = 0

        while len(binary_lst) > 0:
            fixed_list = binary_lst[0]
            binary_lst.remove(fixed_list)
            another_list, is_upper, search_count = self._adjacent_searching_results(fixed_list, binary_lst,
                                                                                    index_bit_length, total_count)
            if search_count >= 0:
                binary_lst.remove(another_list)
                search_counts[search_count] += 1
            else:
                additional += 1
            if is_upper:
                data_set.append(fixed_list)
                data_set.append(another_list)
            else:
                data_set.append(another_list)
                data_set.append(fixed_list)

        results = {}
        for index, count in enumerate(search_counts):
            results[index] = count

        if need_log:
            log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                       "Number of additional bit segment is " + str(additional)
                       + " in original " + str(total_count) + " bit segments.")
            log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                       "In addition, the actual search counts is " + str(results))

        return data_set

    def _adjacent_searching_results(self, fixed_list, other_lists, index_length, total_count, need_check=False):

        random.seed(self.seed)

        if len(other_lists) > 0:
            another_list = other_lists[0]
            n_dna, _ = self._list_to_sequence(fixed_list, another_list)
            return another_list, True, 0

            c_dna, _ = self._list_to_sequence(another_list, fixed_list)
            return another_list, False, 0

        while True:
            # insert at least tenfold interval
            random_index = random.randint(total_count*10, math.pow(2, index_length) - 1)
            index_list = list(map(int, list(str(bin(random_index))[2:].zfill(index_length))))
            print("Random sequence is generating!")

            n_dna, random_list = self._list_to_sequence(fixed_list, index_list)
            return random_list, True, -1

            c_dna, random_list = self._list_to_sequence(index_list, fixed_list)
            return random_list, False, -1

    def _synthesis_sequences(self, data_set, need_log):
        """
        Synthesis sequences by two-dimensional data set.

        :param data_set: original data from file.
        :type: two-dimensional list(int)

        :param need_log: show the log.
        :type: bool
        """

        dna_sequences = []
        for row in range(0, len(data_set), 2):
            if need_log:
                self.monitor.output(row + 2, len(data_set))
            dna_sequence, _ = self._list_to_sequence(data_set[row], data_set[row + 1])
            dna_sequences.append(dna_sequence)

        del data_set

        return dna_sequences

    def _list_to_sequence(self, upper_list, lower_list):
        """
        From two binary sequences to one DNA sequence.

        :param upper_list: the upper binary sequence.
        :type: list

        :param lower_list: the lower binary sequence.
        :type: list
        """

        dna_sequence = []

        if len(upper_list) == len(lower_list):
            for index, (upper_bit, lower_bit) in enumerate(zip(upper_list, lower_list)):
                if upper_bit == 'x' and lower_bit == 'x':
                    dna_sequence.append('T')
                    continue
                if upper_bit == 'x' or lower_bit == 'x':
                    raise Exception \
                        ("The adjacent binary string have different number of next index, which can not be encoded!")

                if index > self.support_spacing:
                    support_base = dna_sequence[index - (self.support_spacing + 1)]
                else:
                    support_base = self.support_bases[index]

                dna_sequence.append(self._binary_to_base(upper_bit, lower_bit, support_base))

            return dna_sequence, None

        addition_length = abs(len(upper_list) - len(lower_list))

        if len(upper_list) > len(lower_list):
            print(len(upper_list), len(lower_list))
            flag = -1
            re_upper_list = copy.deepcopy(upper_list)
            re_lower_list = copy.deepcopy(lower_list + [-1 for _ in range(addition_length)])
        else:
            print(len(upper_list), len(lower_list))
            flag = 1
            re_upper_list = copy.deepcopy(upper_list + [-1 for _ in range(addition_length)])
            re_lower_list = copy.deepcopy(lower_list)

        for index, (upper_bit, lower_bit) in enumerate(zip(re_upper_list, re_lower_list)):
            if upper_bit == 'x' and lower_bit == 'x':
                dna_sequence.append('T')
                continue
            if upper_bit == 'x' or lower_bit == 'x':
                raise Exception \
                    ("The adjacent binary string have different number of next index, which can not be encode!")

            if index > self.support_spacing:
                support_base = dna_sequence[index - (self.support_spacing + 1)]
            else:
                support_base = self.support_bases[index]

            if upper_bit != -1 and lower_bit != -1:
                dna_sequence.append(self._binary_to_base(upper_bit, lower_bit, support_base))
            elif upper_bit == -1:
                is_chosen = False
                for chosen_bit in [0, 1]:
                    current_base = self._binary_to_base(chosen_bit, lower_bit, support_base)
                    re_upper_list[index] = chosen_bit
                    dna_sequence.append(current_base)
                    is_chosen = True
                    break

                if not is_chosen:
                    return None, re_upper_list
            else:
                is_chosen = False
                for chosen_bit in [0, 1]:
                    current_base = self._binary_to_base(upper_bit, chosen_bit, support_base)
                    re_lower_list[index] = chosen_bit
                    dna_sequence.append(current_base)
                    is_chosen = True
                    break

                if not is_chosen:
                    return None, re_lower_list

        if flag == 1:
            return dna_sequence, re_upper_list
        else:
            return dna_sequence, re_lower_list

    def _binary_to_base(self, upper_bit, lower_bit, support_base):
        """
        Get one base from two binary, based on the rules of YYC.

        :param upper_bit: the upper bit, used to identify two of the four bases by RULE 1.
        :type: string

        :param lower_bit: the lower bit, used to identify one of the two bases by RULE 2, after RULE 1.
        :type: string

        :param support_base: the base for support to get base in current position,
                             used to identify one of the two bases by RULE 2, after RULE 1.
        :type: string
        """
        current_options = []
        for index in range(len(self.base_reference)):
            if self.base_reference[index] == int(upper_bit):
                current_options.append(index)

        if self.current_code_matrix[base_index[support_base]][current_options[0]] == int(lower_bit):
            one_base = index_base[current_options[0]]
        else:
            one_base = index_base[current_options[1]]

        return one_base

    # ================================================= decode part ====================================================

    def decode(self, dna_sequences, need_log=False):
        """
        Decode DNA sequences to the data of binary file.

        :param dna_sequences: the DNA sequence of len(matrix) rows.
        :type: list

        :param need_log: show the log.
        :type: bool
        """

        if not dna_sequences:
            log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                       "DNA sequence string set is not existing")

        self.monitor.restore()

        if need_log:
            log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                       "Convert DNA sequences to binary matrix.")

        matrix = self._convert_binaries(dna_sequences, need_log)

        self.monitor.restore()

        return matrix, self.file_size

    def _convert_binaries(self, dna_sequences, need_log):
        """
        Convert DNA sequences to binary matrix.

        :param dna_sequences: the DNA sequence of len(matrix) rows.
        :type: list

        :param need_log: show the log.
        :type: bool
        """

        matrix = []

        last_valid_length = 0
        for row in range(len(dna_sequences)):
            if need_log:
                self.monitor.output(row + 1, len(dna_sequences))

            upper_row_datas, lower_row_datas = self._sequence_to_list(dna_sequences[row])
            matrix.append(upper_row_datas)

            if upper_row_datas != lower_row_datas:
                matrix.append(lower_row_datas)

        del dna_sequences

        return matrix

    def _sequence_to_list(self, dna_sequence, first_idx_length=20, next_idx_length=14):
        """
        Convert one DNA sequence to two-line binary list.

        :param dna_sequence: the DNA sequence of len(matrix) rows.
        :type: list

        :param first_idx_length: the length of the prime index.
        :type: int

        :param next_idx_length: the total length of the minor index.
                                Note that the length of the minor index is the length of the marking plus
                                the valid length.
        :type: int
        """

        upper_row_list = ''
        lower_row_list = ''
        col = 0
        while col < len(dna_sequence):

            col_num = int((col - first_idx_length) / next_idx_length)
            if dna_sequence[col:col+4] == 'TTTT' \
                    and (col - first_idx_length) % next_idx_length == 0 \
                    and [dna_sequence[(first_idx_length + x * next_idx_length):
                    (first_idx_length + x * next_idx_length) + 4] == 'TTTT' for x in range(col_num)] == [1]*col_num:
                upper_row_list += 'xxxx'
                lower_row_list += 'xxxx'
                col += 4

            else:
                if col > self.support_spacing:
                    upper_binary, lower_binary = self._base_to_binary(dna_sequence[col],
                                                                      dna_sequence[col - (self.support_spacing + 1)])
                    upper_row_list += str(upper_binary)
                    lower_row_list += str(lower_binary)
                else:
                    upper_binary, lower_binary = self._base_to_binary(dna_sequence[col], self.support_bases[col])
                    upper_row_list += str(upper_binary)
                    lower_row_list += str(lower_binary)
                col += 1

        return upper_row_list, lower_row_list

    def _base_to_binary(self, current_base, support_base):
        """
        Get two bit from current base and support base, based on the rules of YYC.

        :param current_base: the upper bit, used to identify the upper bit by RULE 1.
        :type: string

        :param support_base: the base for support to get base in current position,
                             used to identify the lower bit by RULE 2.
        :type: string
        """
        upper_bit = self.base_reference[base_index[current_base]]
        lower_bit = self.current_code_matrix[base_index[support_base]][base_index[current_base]]

        return upper_bit, lower_bit
