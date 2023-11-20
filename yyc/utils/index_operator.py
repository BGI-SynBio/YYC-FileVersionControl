import copy
import copyreg
import sys
import numpy as np

import yyc.utils.log as log
from yyc.utils.monitor import Monitor


# noinspection PyProtectedMember
def connect_all(matrix, need_log=False):
    """
    Integrate index and data from the two-dimensional matrix.

    :param matrix: data from input.
    :type: list

    :param need_log: show the log.
    :type: bool
    """
    m = Monitor()
    # index_binary_length = int(len(str(bin(len(matrix)))) - 2)
    index_binary_length = 20

    if need_log:
        log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                   "Add index in the binary matrix.")

    new_matrix = []
    for row in range(len(matrix)):
        if need_log:
            m.output(row + 1, len(matrix))
        new_matrix.append(connect(row, matrix[row], index_binary_length))

    m.restore()

    del matrix, m

    return new_matrix


def connect(index, data, index_binary_length):
    """
    Integrate index and data, list 0100+111101010.

    :param index: the index of data.
    :type: int

    :param data: data from input.
    :type: list

    :param index_binary_length: length of binary index.
    :type: int
    """
    bin_index = list(map(int, list(str(bin(index))[2:].zfill(index_binary_length))))
    one_list = bin_index + data

    return one_list


def divide(one_list, index_binary_length):
    """
    Separate data from the index in a binary sequence.

    :param one_list: one binary sequence.
    :type: list

    :param index_binary_length: length of binary index.
    :type: int
    """
    # Convert binary index to decimal.
    index = int("".join(list(map(str, one_list[:index_binary_length]))), 2)
    data = one_list[index_binary_length:]

    return index, data


def divide_all_version(binary_lst, need_log=False, first_idx_length=20, next_idx_length=14):
    """
    Separate data from indexes in binary sequences.

    :param binary_lst: the binary sequences.
    :type: list

    :param need_log: show the log.
    :type: bool

    :param first_idx_length: the length of the prime index.
    :type: int

    :param next_idx_length: the total length of the minor index.
                            Note that the length of the minor index is the length of the marking plus the valid length.
    :type: int
    """
    m = Monitor()
    # index_binary_length = int(len(str(bin(len(matrix)))) - 2)
    # index_binary_length = 20

    if need_log:
        log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                   "Divide index and data from binary matrix.")

    indexs = []
    datas = []
    indexs_int = []
    datas_int = []
    have_next_index_count = 0
    for idx, binary_string in enumerate(binary_lst):
        if need_log:
            m.output(idx + 1, len(binary_lst))
        if 'xxxx' not in binary_string:
            index = binary_string[:first_idx_length]
            data = binary_string[first_idx_length:]
        else:
            next_idx_num = binary_string.count('xxxx')
            index = binary_string[:first_idx_length+next_idx_num*next_idx_length]
            data = binary_string[first_idx_length+next_idx_num*next_idx_length:]
            have_next_index_count += 1

        indexs.append(index)
        datas.append(data)
        indexs_int.append([int(x, 2) for x in index.split('xxxx')])

    return indexs, datas, indexs_int


# noinspection PyProtectedMember
def sort_order_version(indexs, datas, indexs_binary, need_log=False):
    """
    Restore data in order of index.

    :param indexes: the indexes of datas. For example, 20-1-2.
    :type: list

    :param datas: the disordered datas, the locations of these are corresponding to parameter "index".

    :param indexs_binary: the binary indexes of datas.
    :type: list

    :param need_log: show the log.
    :type: bool
    """
    m = Monitor()

    if need_log:
        log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                   "Restore data order according to index.")

    max_idx_num = max(len(x) for x in indexs)
    fill_indexs = list(map(lambda x:x+[0]*(max_idx_num-len(x)), indexs))
    sort_indexs = sorted(fill_indexs, key=lambda x: ([x[i] for i in range(max_idx_num)]))
    sort_datas = []
    sort_indexs_datas = []
    for index in sort_indexs:
        idx = fill_indexs.index(index)
        sort_datas.append(datas[idx])
        sort_indexs_datas.append(indexs_binary[idx] + datas[idx])

    m.restore()
    del indexs, datas, indexs_binary, m

    return sort_datas, sort_indexs_datas


def multiple_index_sort(max_index_num, min_index_num, multiple_index):

    last_index_items = list(multiple_index.items())
    num = min_index_num + 1

    while num <= max_index_num:
        need_sorted_idx = []
        for idx, temp in enumerate(last_index_items):
            if len(temp[1]) > num:
                need_sorted_idx.append(idx)

        sorted_idx = copy.deepcopy(need_sorted_idx)
        sorted_index_items = copy.deepcopy(last_index_items)
        sorted_idx = sorted(need_sorted_idx, key=lambda x:([x[1][i] for i in range(num)]))
        for idx in need_sorted_idx:
            sorted_index_items[idx] = last_index_items[sorted_idx[idx]]

        num += 1
        last_index_items = copy.deepcopy(sorted_index_items)

    return last_index_items
