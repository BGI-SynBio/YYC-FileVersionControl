import struct
import math
import sys
import os

import yyc.utils.log as log
from yyc.utils.monitor import Monitor


# noinspection PyProtectedMember
def read_binary_from_all(path, payload_length, need_log=False):
    """
    Reading binary matrix from document.

    :param path: file path.
    :type: string

    :param payload_length: the binary segment length used for DNA sequence generation.
                           Considering current DNA synthesis technique limitation,
                           we usually set 120 as default segment length.
    :type: int

    :param need_log: show the log.
    :type: bool
    """

    m = Monitor()
    try:

        # Open selected file
        with open(path, mode="rb") as file:

            if need_log:
                log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                           "Read binary matrix from file: " + path)

            size = os.path.getsize(path)

            # Set init storage matrix
            matrix = [[0 for _ in range(payload_length)] for _ in range(math.ceil(size * 8 / payload_length))]

            row = 0
            col = 0
            for byte_index in range(size):
                if need_log:
                    m.output(byte_index + 1, size)
                # Read a file as bytes
                one_byte = file.read(1)
                element = list(map(int, list(str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8))))
                for bit_index in range(8):
                    matrix[row][col] = element[bit_index]
                    col += 1
                    if col == payload_length:
                        col = 0
                        row += 1

        if int(len(str(bin(len(matrix)))) - 2) * 7 > payload_length:
            if need_log:
                log.output(log.WARN, str(__name__), str(sys._getframe().f_code.co_name),
                           "The proportion of index in whole sequence may be high. \n"
                           "It is recommended to increase the length of output DNA sequences "
                           "or to divide the file into more segment pools")

        return matrix, size
    except IOError:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The file selection operation was not performed correctly. Please execute the operation again!")


def read_binary_string_from_all(path, need_log=False):

    m = Monitor()
    try:

        # Open selected file
        with open(path, mode="rb") as file:

            if need_log:
                log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                           "Read binary matrix from file: " + path)

            size = os.path.getsize(path)
            binary_string = ""

            for byte_index in range(size):
                if need_log:
                    m.output(byte_index + 1, size)
                # Read a file as bytes
                one_byte = file.read(1)
                element = str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8)
                binary_string += element

        return binary_string, size

    except IOError:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The file selection operation was not performed correctly. Please execute the operation again!")


# noinspection PyBroadException,PyProtectedMember
def write_all_from_binary(path, matrix, size, need_log=False):
    """
    Writing binary matrix to document.

    :param path: file path.
    :type: string

    :param matrix: a matrix in which each row represents a binary sequence that will be used for DNA sequence generation.
    :type: list

    :param size: this refers to file size, to reduce redundant bits when transferring DNA to binary files.
    :type: int

    :param need_log: show the log.
    :type: bool
    """
    m = Monitor()

    try:
        with open(path, "wb+") as file:
            if need_log:
                log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                           "Write file from binary matrix: " + path)

            bit_index = 0
            temp_byte = 0
            count = 0
            for row in range(len(matrix)):
                print(matrix[row])
                if need_log:
                    m.output(row + 1, len(matrix))

                col = 0
                while col < len(matrix[row]):
                    if col == 0 and matrix[row][:15] == '0'*15:
                        count += 1
                        need_extract = int('0b'+matrix[row][15:22], 2)
                        for c in range(len(matrix[row])-need_extract, len(matrix[row])):
                            bit_index += 1
                            temp_byte *= 2
                            temp_byte += int(matrix[row][c])
                            if bit_index == 8:
                                file.write(struct.pack("B", int(temp_byte)))
                                bit_index = 0
                                temp_byte = 0
                        col = len(matrix[row])
                    else:
                        bit_index += 1
                        temp_byte *= 2
                        temp_byte += int(matrix[row][col])
                        if bit_index == 8:
                            file.write(struct.pack("B", int(temp_byte)))
                            bit_index = 0
                            temp_byte = 0
                        col += 1

    except IOError:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The file selection operation was not performed correctly. Please execute the operation again!")


# noinspection PyBroadException,PyProtectedMember
def read_dna_file(path, need_log=False):
    """
    Reading DNA sequence set from documents.

    :param path: file path.
    :type: string

    :param need_log: show the log.
    :type: bool
    """

    m = Monitor()

    dna_sequences = []
    dna_sequences2 = []

    try:
        with open(path, "r") as file:
            if need_log:
                log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                           "Read DNA sequences from file: " + path)

            # Read current file by line
            lines = file.readlines()
            for index in range(len(lines)):
                if need_log:
                    m.output(index + 1, len(lines))
                line = lines[index]
                dna_sequences.append([line[col] for col in range(len(line) - 1)])
                dna_sequences2.append(''.join([line[col] for col in range(len(line) - 1)]))

        return dna_sequences, dna_sequences2
    except IOError:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The file selection operation was not performed correctly. Please execute the operation again!")


# noinspection PyProtectedMember,PyBroadException
def write_dna_file(path, dna_sequences, need_log=False):
    """
    Writing DNA sequence set to documents.

    :param path: file path.
    :type: string

    :param dna_sequences: generated DNA sequences.
    :type: one-dimensional list(string)

    :param need_log: show the log.
    :type: bool
    """

    m = Monitor()

    try:
        with open(path, "w") as file:
            if need_log:
                log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                           "Write DNA sequences to file: " + path)
            for row in range(len(dna_sequences)):
                if need_log:
                    m.output(row + 1, len(dna_sequences))
                file.write("".join(dna_sequences[row]) + "\n")
        return dna_sequences
    except IOError:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The file selection operation was not performed correctly. Please execute the operation again!")
