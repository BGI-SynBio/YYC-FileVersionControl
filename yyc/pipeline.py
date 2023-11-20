import copy
import sys
import math
from yyc.utils import log, data_handle, index_operator, model_saver
import numpy as np


# noinspection PyProtectedMember
def encode_original(method, input_path, output_path,
                    model_path=None, verify=None, need_index=True, payload_length=120, need_log=False):
    """
    Use the selected method, convert the original version file to DNA sequences and output the DNA sequences to a file.

    :param method: transcoding method.
                   If you have model file, you can use this function without method.
    :type: object

    :param input_path: the path of the original version file.
    :type: string

    :param output_path: the path of the output file.
    :type: string

    :param model_path: the path of model file if you want to save.
    :type: string

    :param verify: error correction method.
    :type: object

    :param need_index: whether need index for the DNA sequences.
    :type: bool

    :param payload_length: the binary segment length used for DNA sequence generation.
                           Considering current DNA synthesis technique limitation,
                           we usually set 120 as default segment length.
    :type: int

    :param need_log: show the log.
    :type: bool
    """

    if input_path is None or len(input_path) == 0:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The input file path is invalid!")

    if output_path is None or len(input_path) == 0:
        log.output(log.ERROR, str(__name__), str(sys._getframe().f_code.co_name),
                   "The output file path is invalid!")

    input_matrix, size = data_handle.read_binary_from_all(input_path, payload_length, need_log)

    if need_index:
        input_matrix = index_operator.connect_all(input_matrix, need_log)

    if verify is not None:
        input_matrix = verify.add_for_matrix(input_matrix, need_log)

    dna_sequences = method.encode(input_matrix, need_log)

    if model_path is not None:
        model_saver.save_model(model_path, {"method": method, "verify": verify})

    data_handle.write_dna_file(output_path, dna_sequences, need_log)


# noinspection PyProtectedMember
def encode(method, modifications, modified_files, match_files, last_version_dnas_path, output_path, need_index=True,
           segment_length=140, first_idx_length=20, next_idx_length=14, zero_mark=22, limit_next_index_num=6,
           need_log=False):
    """
    Use the selected method, encode the current version file based on the dna file of the last version and output the
    dna file of current version.

    :param method: transcoding method.
                   If you have model file, you can use this function without method.
    :type: object

    :param modifications: type of the modification operation, including insert, delete and replace.
    :type: list(char)

    :param modified_files: content of the modification operation.
                           Note that replace operation is considered as a delete operation plus a insert operation,
                           so the modified_files of the replace operation contains two files.
    :type: list(char or list(char))

    :param match_files: the content in front of the modified part in order to locate the modified part.
                        Note that it means the modified part is in the beginning of the text if there is no match file.
    :type: list(char)

    :param last_version_dnas_path: the path of the last version dna files.
    :type: string

    :param output_path: the path of the output file.
    :type: string

    :param need_index: whether need index for the DNA sequences.
    :type: bool

    :param segment_length:  The cut length of DNA sequence plus the length of the prime index.
    :type: int

    :param first_idx_length: the length of the prime index.
    :type: int

    :param next_idx_length: the total length of the minor index.
                            Note that the length of the minor index is the length of the marking plus the valid length.
    :type: int

    :param zero_mark: the number of zero for marking the binary sequences that need to add invalid bianry stream.
    :type: int

    :param limit_next_index_num: the valid length of the minor index.
    :type: int

    :param need_log: show the log.
    :type: bool
    """

    _, last_version_dnas = data_handle.read_dna_file(last_version_dnas_path, need_log)
    last_version_matrix, size = method.decode(last_version_dnas, need_log)
    # print(last_version_matrix)

    # identify the 'xxxx' to minor index, such as "20-1-1"
    last_version_binary_indexes, last_version_data_set, last_version_indexes = \
        index_operator.divide_all_version(last_version_matrix, need_log)
    new_version_dnas = copy.deepcopy(last_version_matrix)

    last_version_idx_len = {}
    last_version_binary_string = ""
    valid_last_version_data_set = []
    cur_sum_count = 0
    for idx, data in enumerate(last_version_data_set):
        if data[:15] == '0'*15:
            print(data[:15])
            need_extract = int('0b' + data[15:zero_mark], 2)
            valid_data = data[-need_extract:]
        else:
            valid_data = data

        cur_sum_count += len(valid_data)
        valid_last_version_data_set.append(valid_data)
        last_version_idx_len[idx] = [cur_sum_count, len(valid_data)]
        last_version_binary_string += "".join([str(d) for d in valid_data])

    last_version_data_set = copy.deepcopy(valid_last_version_data_set)
    print(modifications)

    for idx, modification in enumerate(modifications):
        modified_file = modified_files[idx]
        print(modified_file)
        if modification == 'replace':
            modified_binary_string = {'delete': '', 'insert': ''}
            modified_binary_string['delete'], _ = data_handle.read_binary_string_from_all(modified_file[0], need_log)
            modified_binary_string['insert'], _ = data_handle.read_binary_string_from_all(modified_file[1], need_log)
        else:
            modified_binary_string, _ = data_handle.read_binary_string_from_all(modified_file, need_log)

        match_file = match_files[idx]
        print(match_file)
        if match_file != "":
            match_binary_string, match_size = data_handle.read_binary_string_from_all(match_file, need_log)
        else:
            match_binary_string = False

        # match to obtain the modified start dna index
        additional_binary = {'head': '', 'end': ''}
        if match_binary_string:
            match_start_str = last_version_binary_string.index(match_binary_string)
            for idx, temp in last_version_idx_len.items():
                str_count = temp[0]
                if str_count < match_start_str + 1:
                    pass
                else:
                    match_start_dna = idx
                    break
            match_end_str = match_start_str + len(match_binary_string) - 1

            # get last_index, next_index, additional_binary
            additional_binary = {'head': '', 'end': ''}
            last_index = False
            next_index = False
            last_binary_index = False
            next_binary_index = False
            for idx, temp in enumerate(list(last_version_idx_len.values())[match_start_dna:]):
                sum_count = temp[0]
                if sum_count < match_end_str + 1:
                    pass
                elif sum_count == match_end_str + 1:
                    match_end_dna = match_start_dna + idx
                    if (match_end_dna + 1) % 2 == 0:
                        if modification == 'insert':
                            if match_end_dna >= 0:
                                last_index = last_version_indexes[match_end_dna]
                                last_binary_index = last_version_binary_indexes[match_end_dna]
                            if match_end_dna < len(last_version_indexes) - 1:
                                next_index = last_version_indexes[match_end_dna + 1]
                                next_binary_index = last_version_binary_indexes[match_end_dna + 1]

                        else:
                            if modification == 'delete':
                                modified_end_str = match_end_str + len(modified_binary_string)
                            else:
                                modified_end_str = match_end_str + len(modified_binary_string['delete'])

                            for id, x in enumerate(list(last_version_idx_len.values())[match_end_dna + 1:]):
                                count = x[0]
                                modified_end_dna = match_end_dna + id + 1
                                new_version_dnas.remove(last_version_matrix[modified_end_dna])
                                if count < modified_end_str + 1:
                                    pass
                                elif count == modified_end_str + 1:
                                    break
                                else:
                                    additional_count = count - (modified_end_str + 1)
                                    additional_string_right = \
                                        last_version_data_set[modified_end_dna][-additional_count:]
                                    additional_binary['end'] = additional_string_right
                                    break

                            if (modified_end_dna - match_end_dna) % 2 == 0:
                                if match_end_dna >= 0:
                                    last_index = last_version_indexes[match_end_dna]
                                    last_binary_index = last_version_binary_indexes[match_end_dna]
                                if modified_end_dna < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 1]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 1]
                                pass
                            else:
                                if modified_end_dna < len(last_version_data_set) - 1:
                                    additional_binary['end'] = additional_binary['end'] \
                                                               + last_version_data_set[modified_end_dna + 1]
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna + 1])
                                if match_end_dna >= 0:
                                    last_index = last_version_indexes[match_end_dna]
                                    last_binary_index = last_version_binary_indexes[match_end_dna]
                                if modified_end_dna < len(last_version_indexes) - 2:
                                    next_index = last_version_indexes[modified_end_dna + 2]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 2]
                    else:
                        additional_binary['head'] = last_version_data_set[match_end_dna]
                        new_version_dnas.remove(last_version_matrix[match_end_dna])
                        if modification == 'insert':
                            if match_end_dna + 1 <= len(last_version_data_set) - 1:
                                additional_binary['end'] = last_version_data_set[match_end_dna + 1]
                                new_version_dnas.remove(last_version_matrix[match_end_dna + 1])

                            if match_end_dna > 0:
                                last_index = last_version_indexes[match_end_dna - 1]
                                last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                            if match_end_dna < len(last_version_indexes) - 2:
                                next_index = last_version_indexes[match_end_dna + 2]
                                next_binary_index = last_version_binary_indexes[match_end_dna + 2]
                        else:
                            if modification == 'delete':
                                modified_end_str = match_end_str + len(modified_binary_string)
                            else:
                                modified_end_str = match_end_str + len(modified_binary_string['delete'])

                            for id, x in enumerate(list(last_version_idx_len.values())[match_end_dna+1:]):
                                count = x[0]
                                modified_end_dna = match_end_dna + id + 1
                                new_version_dnas.remove(last_version_matrix[modified_end_dna])
                                if count < modified_end_str + 1:
                                    pass
                                elif count == modified_end_str + 1:
                                    break
                                else:
                                    additional_count = count - (modified_end_str + 1)
                                    additional_string_right = \
                                        last_version_data_set[modified_end_dna][-additional_count:]
                                    additional_binary['end'] = additional_string_right
                                    # new_version_dnas.remove(last_version_dnas[modified_end_dna])
                                    break

                            if (modified_end_dna - match_end_dna) % 2 == 1:
                                # additional_binary.append(len(additional_binary[0])*'0')
                                # 需要additional_binary最终添加索引后变成偶数
                                if match_end_dna > 0:
                                    last_index = last_version_indexes[match_end_dna - 1]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                                if modified_end_dna < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 1]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 1]
                                pass
                            else:
                                if modified_end_dna < len(last_version_data_set) - 1:
                                    additional_binary['end'] = additional_binary['end'] \
                                                               + last_version_data_set[modified_end_dna + 1]
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna + 1])
                                if match_end_dna > 0:
                                    last_index = last_version_indexes[match_end_dna - 1]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                                if modified_end_dna + 1 < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 2]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 2]
                    break

                else:
                    match_end_dna = match_start_dna + idx
                    additional_count = sum_count - (match_end_str + 1)
                    additional_string_left = last_version_data_set[match_end_dna][:-additional_count]
                    additional_string_right = last_version_data_set[match_end_dna][-additional_count:]
                    if (match_end_dna + 1) % 2 == 0:
                        additional_binary['head'] = last_version_data_set[match_end_dna - 1] + additional_string_left

                        new_version_dnas.remove(last_version_matrix[match_end_dna - 1])
                        new_version_dnas.remove(last_version_matrix[match_end_dna])

                        if modification == 'insert':
                            additional_binary['end'] = additional_string_right

                            if match_end_dna - 2 >= 0:
                                last_index = last_version_indexes[match_end_dna - 2]
                                last_binary_index = last_version_binary_indexes[match_end_dna - 2]
                            if match_end_dna < len(last_version_indexes) - 1:
                                next_index = last_version_indexes[match_end_dna + 1]
                                next_binary_index = last_version_binary_indexes[match_end_dna + 1]

                        else:
                            if modification == 'delete':
                                modified_end_str = match_end_str + len(modified_binary_string)
                            else:
                                modified_end_str = match_end_str + len(modified_binary_string['delete'])

                            for id, x in enumerate(list(last_version_idx_len.values())[match_end_dna:]):
                                count = x[0]
                                modified_end_dna = match_end_dna + id

                                if id > 0:
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna])
                                if count < modified_end_str + 1:
                                    pass
                                elif count == modified_end_str + 1:
                                    break
                                else:
                                    additional_count = count - (modified_end_str + 1)
                                    additional_string_right = \
                                        last_version_data_set[modified_end_dna][-additional_count:]
                                    additional_binary['end'] = additional_string_right
                                    break
                            if (modified_end_dna - match_end_dna) % 2 == 1:
                                if modified_end_dna < len(last_version_data_set) - 1:
                                    additional_binary['end'] = additional_binary['end'] \
                                                               + last_version_data_set[modified_end_dna + 1]
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna + 1])
                                if match_end_dna - 2 >= 0:
                                    last_index = last_version_indexes[match_end_dna - 2]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 2]
                                if modified_end_dna + 1 < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 2]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 2]
                                pass
                            else:
                                if match_end_dna - 2 >= 0:
                                    last_index = last_version_indexes[match_end_dna - 2]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 2]
                                if modified_end_dna < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 1]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 1]
                    else:
                        additional_binary['head'] = additional_string_left
                        new_version_dnas.remove(last_version_matrix[match_end_dna])

                        if modification == 'insert':
                            additional_binary['end'] = additional_string_right
                            if match_end_dna + 1 <= len(last_version_data_set) - 1:
                                additional_binary['end'] = additional_binary['end'] \
                                                           + last_version_data_set[match_end_dna + 1]
                                new_version_dnas.remove(last_version_matrix[match_end_dna + 1])
                            if match_end_dna - 1 >= 0:
                                last_index = last_version_indexes[match_end_dna - 1]
                                last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                            if match_end_dna + 1 < len(last_version_indexes) - 1:
                                next_index = last_version_indexes[match_end_dna + 2]
                                next_binary_index = last_version_binary_indexes[match_end_dna + 2]

                        else:
                            if modification == 'delete':
                                modified_end_str = match_end_str + len(modified_binary_string)
                            else:
                                modified_end_str = match_end_str + len(modified_binary_string['delete'])

                            for id, x in enumerate(list(last_version_idx_len.values())[match_end_dna:]):
                                count = x[0]
                                modified_end_dna = match_end_dna + id
                                if id > 0:
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna])
                                if count < modified_end_str + 1:
                                    pass
                                elif count == modified_end_str + 1:
                                    break
                                else:
                                    additional_count = count - (modified_end_str + 1)
                                    additional_string_right = \
                                        last_version_data_set[modified_end_dna][-additional_count:]
                                    additional_binary['end'] = additional_string_right
                                    break
                            if (modified_end_dna - match_end_dna) % 2 == 1:
                                if match_end_dna > 0:
                                    last_index = last_version_indexes[match_end_dna - 1]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                                if modified_end_dna < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 1]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 1]
                                pass
                            else:
                                if modified_end_dna < len(last_version_data_set) - 1:
                                    additional_binary['end'] = additional_binary['end'] \
                                                               + last_version_data_set[modified_end_dna + 1]
                                    new_version_dnas.remove(last_version_matrix[modified_end_dna + 1])
                                if match_end_dna > 0:
                                    last_index = last_version_indexes[match_end_dna - 1]
                                    last_binary_index = last_version_binary_indexes[match_end_dna - 1]
                                if modified_end_dna + 1 < len(last_version_indexes) - 1:
                                    next_index = last_version_indexes[modified_end_dna + 2]
                                    next_binary_index = last_version_binary_indexes[modified_end_dna + 2]
                    break

        else:
            last_index = False
            last_binary_index = False
            next_index = last_version_indexes[0]
            next_binary_index = last_version_binary_indexes[0]

        if modification == 'replace':
            binary_string = additional_binary['head'] + modified_binary_string['insert'] + additional_binary['end']
        elif modification == 'insert':
            binary_string = additional_binary['head'] + modified_binary_string + additional_binary['end']
        else:
            binary_string = additional_binary['head'] + additional_binary['end']

        add_end_length = len(additional_binary['end'])
        end_remainder = add_end_length % 8
        modified_binary_lst = []
        if need_index:
            limited_idx_num = pow(2, next_idx_length - 4) - 1
            if last_index and next_index:
                if len(last_index) > len(next_index):
                    first_unequal = len(next_index)
                    for i in range(len(next_index)):
                        if last_index[i] != next_index[i]:
                            first_unequal = i
                            break
                else:
                    first_unequal = len(last_index)
                    for i in range(len(last_index)):
                        if last_index[i] != next_index[i]:
                            first_unequal = i
                            break
                internal = len(last_index) - first_unequal - 1
                if internal == 0:
                    need_zero = 0
                    add_next_idx = [last_index[-1]+1]
                    next_idx_num = len(last_index) - 1
                    payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                    modified_count = math.ceil(len(binary_string) / payload_length)
                    last_next_diff = next_index[len(last_index) - 1] - last_index[len(last_index) - 1]
                    if modified_count+last_index[-1]+3 >= last_next_diff:
                        need_zero = 1
                        add_next_idx = [last_index[-1], 1]
                        next_idx_num = len(last_index) - 1 + 1
                        payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                        modified_count = math.ceil(len(binary_string) / payload_length)
                        if modified_count+3 > limited_idx_num * last_next_diff:
                            need_zero = 2
                            add_next_idx = [1, 1]
                            next_idx_num = len(last_index) - 1 + 2
                            payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                            modified_count = math.ceil(len(binary_string) / payload_length)
                            if math.ceil((modified_count+3) / limited_idx_num) > limited_idx_num:
                                raise Exception("Need more number of next index!")

                elif internal == 1:
                    need_zero = 0
                    add_next_idx = [last_index[-2], last_index[-1]+1]
                    next_idx_num = len(last_index) - 1
                    payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                    modified_count = math.ceil(len(binary_string) / payload_length)
                    last_next_diff = next_index[len(last_index)-2] - last_index[-2]
                    if modified_count + last_index[-1] + 3 > limited_idx_num * last_next_diff:
                        need_zero = 1
                        add_next_idx = [1, last_index[-1]+1]
                        next_idx_num = len(last_index) - 1 + 1
                        payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                        modified_count = math.ceil(len(binary_string) / payload_length)
                        if math.ceil((modified_count + last_index[-1] - 1 + 3) / limited_idx_num) > limited_idx_num:
                            raise Exception("Need more number of next index!")

                else:
                    need_zero = 0
                    add_next_idx = [last_index[-2], last_index[-1]+1]
                    next_idx_num = len(last_index) - 1
                    payload_length = segment_length - (next_idx_num * next_idx_length + first_idx_length)
                    modified_count = math.ceil(len(binary_string) / payload_length)
                    if math.ceil((modified_count + last_index[-1] + 3) / limited_idx_num) > limited_idx_num - last_index[-2] + 1:
                        raise Exception("Need more number of next index!")

                if next_idx_num > limit_next_index_num:
                    raise Exception("The number of the next index exceed its limitation!")
                last_string_length = len(binary_string) % payload_length

                if last_string_length > 0:
                    if payload_length - last_string_length >= zero_mark:
                        mark_num = bin(last_string_length)[2:].zfill(7)
                        binary_string = binary_string[:-last_string_length] + 15 * '0' + mark_num \
                                        + (payload_length - last_string_length - zero_mark) * '0' \
                                        + binary_string[-last_string_length:]
                        if payload_length - last_string_length > zero_mark:
                            pass

                    else:
                        mark_num = bin(payload_length - zero_mark)[2:].zfill(7)
                        mark_num_int = int('0b'+mark_num, 2)

                        binary_string_adding = binary_string[:-last_string_length] + 15 * '0' + mark_num \
                                               + binary_string[-last_string_length: -last_string_length+mark_num_int]
                        remain_length = last_string_length - mark_num_int
                        while remain_length:
                            modified_count += 1
                            if remain_length <= mark_num_int:
                                mark_num = bin(remain_length)[2:].zfill(7)
                                binary_string_adding = binary_string_adding + 15 * '0' \
                                                       + mark_num + (mark_num_int - remain_length) * '0' \
                                                       + binary_string[-remain_length:]
                                remain_length = 0
                            else:
                                binary_string_adding = binary_string_adding + 15 * '0' + mark_num \
                                                       + binary_string[-remain_length:-remain_length+mark_num_int]
                                remain_length = remain_length - mark_num_int
                        binary_string = copy.deepcopy(binary_string_adding)

                if modified_count % 2 == 0:
                    pass
                else:
                    binary_string = binary_string + '0' * payload_length
                    modified_count += 1

                for c in range(modified_count):
                    if next_idx_num == len(last_index) - 1:
                        if len(add_next_idx) == 1:
                            if len(last_index) == 1:
                                binary_sequence = bin(add_next_idx[0]+c)[2:].zfill(first_idx_length) \
                                                  + binary_string[c*payload_length:(c+1)*payload_length]
                            else:
                                binary_sequence = last_binary_index[:-next_idx_length] + 'xxxx' \
                                                  + bin(add_next_idx[0]+c)[2:].zfill(next_idx_length - 4) \
                                                  + binary_string[c*payload_length:(c+1)*payload_length]
                        elif len(add_next_idx) == 2:
                            num = math.ceil((c + add_next_idx[-1]) / limited_idx_num) - 1
                            final_idx = c + add_next_idx[-1] - limited_idx_num * num

                            if len(last_index) == 2:
                                binary_sequence = bin(add_next_idx[0] + num)[2:].zfill(first_idx_length) \
                                                  + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                                  + binary_string[c * payload_length:(c + 1) * payload_length]
                            else:
                                binary_sequence = last_binary_index[:-(next_idx_length*2)] + 'xxxx' \
                                                  + bin(add_next_idx[0] + num)[2:].zfill(next_idx_length - 4) \
                                                  + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                                  + binary_string[c * payload_length:(c + 1) * payload_length]

                    elif next_idx_num == len(last_index):
                        if len(add_next_idx) == 1:
                            binary_sequence = last_binary_index + 'xxxx' \
                                              + bin(add_next_idx[0] + c)[2:].zfill(next_idx_length - 4) \
                                              + binary_string[c * payload_length:(c + 1) * payload_length]
                        elif len(add_next_idx) == 2:
                            num = math.ceil((c + add_next_idx[-1]) / limited_idx_num) - 1
                            final_idx = c + add_next_idx[-1] - limited_idx_num * num
                            if len(last_index) == 1:
                                binary_sequence = bin(add_next_idx[0] + num)[2:].zfill(first_idx_length) \
                                                  + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                                  + binary_string[c * payload_length:(c + 1) * payload_length]
                                print(add_next_idx[0] + num, final_idx)
                            else:
                                binary_sequence = last_binary_index[:-next_idx_length] + 'xxxx' \
                                                  + bin(add_next_idx[0] + num)[2:].zfill(next_idx_length - 4) \
                                                  + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                                  + binary_string[c * payload_length:(c + 1) * payload_length]

                    elif next_idx_num == len(last_index) + 1:
                        num = math.ceil((c + add_next_idx[-1]) / limited_idx_num) - 1
                        final_idx = c + add_next_idx[-1] - limited_idx_num * num
                        binary_sequence = last_binary_index + 'xxxx' \
                                          + bin(add_next_idx[0] + num)[2:].zfill(next_idx_length - 4) \
                                          + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                          + binary_string[c * payload_length:(c + 1) * payload_length]
                    modified_binary_lst.append(binary_sequence)

            elif last_index:
                payload_length = segment_length - first_idx_length
                modified_count = math.ceil(len(binary_string) / payload_length)
                if modified_count + 3 <= pow(2, first_idx_length-4) - 1 - last_index[0]:
                    next_idx_num = 0
                    add_next_idx = [last_index[0]+1]
                else:
                    next_idx_num = 1
                    payload_length = segment_length - first_idx_length - next_idx_length*next_idx_num
                    modified_count = math.ceil(len(binary_string) / payload_length)
                    if len(last_index) == 1:
                        add_next_idx = [last_index[0], 1]
                    else:
                        add_next_idx = [last_index[0], last_index[1]+1]
                    if (modified_count + add_next_idx[-1] - 1 + 3) / limited_idx_num > \
                            pow(2, first_idx_length) - 1 - add_next_idx[0] + 1:
                        raise Exception("Need more number of next index!")

                if len(next_idx_num) > limit_next_index_num:
                    raise Exception("The number of the next index exceed its limitation!")

                last_string_length = len(binary_string) % payload_length

                if last_string_length > 0:
                    if payload_length - last_string_length >= zero_mark:
                        binary_string = binary_string[:-last_string_length] + binary_string[-last_string_length:] \
                                        + (payload_length - last_string_length) * '0'
                        if payload_length - last_string_length > zero_mark:
                            pass

                    else:
                        mark_num = bin(payload_length - zero_mark)[2:].zfill(7)
                        mark_num_int = int('0b' + mark_num, 2)
                        binary_string_adding = binary_string[:-last_string_length] + 15 * '0' + mark_num \
                                               + binary_string[-last_string_length: -last_string_length + mark_num_int]
                        remain_length = last_string_length - mark_num_int

                        while remain_length:
                            modified_count += 1
                            if remain_length <= mark_num_int:
                                mark_num = bin(remain_length)[2:].zfill(7)
                                binary_string_adding = binary_string_adding + 15 * '0' \
                                                       + mark_num + (mark_num_int - remain_length) * '0' \
                                                       + binary_string[-remain_length:]
                                remain_length = 0
                            else:
                                binary_string_adding = binary_string_adding + 15 * '0' + mark_num \
                                                       + binary_string[-remain_length:-remain_length + mark_num_int]
                                remain_length = remain_length - mark_num_int
                        binary_string = copy.deepcopy(binary_string_adding)

                if modified_count % 2 == 0:
                    pass
                else:
                    binary_string = binary_string + '0' * payload_length
                    modified_count += 1

                for c in range(modified_count):
                    if next_idx_num == 0:
                        binary_sequence = bin(c + add_next_idx[0])[2:].zfill(first_idx_length) \
                                          + binary_string[c * payload_length:(c + 1) * payload_length]
                        modified_binary_lst.append(binary_sequence)
                    elif next_idx_num == 1:
                        num = math.ceil((c + add_next_idx[-1]) / limited_idx_num) - 1
                        final_idx = c + add_next_idx[-1] - limited_idx_num * num
                        binary_sequence = bin(add_next_idx[0] + num)[2:].zfill(first_idx_length) \
                                          + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                          + binary_string[c * payload_length:(c + 1) * payload_length]
                        modified_binary_lst.append(binary_sequence)

            else:
                need_zero = False
                for idx, s in enumerate(next_index):
                    s = int(s)
                    if s > 0:
                        need_zero = idx
                        break

                next_idx_num = need_zero + 1
                add_next_idx = [1]
                payload_length = segment_length - first_idx_length - next_idx_length * next_idx_num
                modified_count = math.ceil(len(binary_string) / payload_length)
                if modified_count + 3 > limited_idx_num:
                    need_zero += 1
                    next_idx_num = need_zero + 1
                    add_next_idx = [1, 1]
                    payload_length = segment_length - first_idx_length - next_idx_length * next_idx_num
                    modified_count = math.ceil(len(binary_string) / payload_length)
                    if (modified_count + 3) / limited_idx_num > limited_idx_num:
                        raise Exception("Need more number of next index!")

                if next_idx_num > limit_next_index_num:
                    raise Exception("The number of the next index exceed its limitation!")
                last_string_length = len(binary_string) % payload_length

                if last_string_length > 0:
                    if payload_length - last_string_length >= zero_mark:
                        binary_string = binary_string[:-last_string_length] + binary_string[-last_string_length:] \
                                        + (payload_length - last_string_length) * '0'
                        if payload_length - last_string_length > zero_mark:
                            pass

                    else:
                        mark_num = bin(payload_length - zero_mark)[2:].zfill(7)
                        mark_num_int = int('0b' + mark_num, 2)

                        binary_string_adding = binary_string[:-last_string_length] + 15 * '0' + mark_num \
                                               + binary_string[-last_string_length: -last_string_length + mark_num_int]
                        remain_length = last_string_length - mark_num_int
                        while remain_length:
                            modified_count += 1
                            if remain_length <= mark_num_int:
                                mark_num = bin(remain_length)[2:].zfill(7)
                                binary_string_adding = binary_string_adding + 15 * '0' \
                                                       + mark_num + (mark_num_int - remain_length) * '0' \
                                                       + binary_string[-remain_length:]
                                remain_length = 0
                            else:
                                binary_string_adding = binary_string_adding + 15 * '0' + mark_num \
                                                       + binary_string[-remain_length:-remain_length + mark_num_int]
                                remain_length = remain_length - mark_num_int
                        binary_string = copy.deepcopy(binary_string_adding)

                if modified_count % 2 == 0:
                    pass
                else:
                    binary_string = binary_string + '0' * payload_length
                    modified_count += 1

                for c in range(modified_count):
                    if len(add_next_idx) == 1:
                        binary_sequence = bin(0)[2:].zfill(first_idx_length) + 'xxxx'\
                                          + bin(c+1)[2:].zfill(next_idx_length - 4) \
                                          + binary_string[c * payload_length:(c + 1) * payload_length]
                        modified_binary_lst.append(binary_sequence)
                    else:
                        num = math.ceil((c + add_next_idx[-1]) / limited_idx_num) - 1
                        final_idx = c + add_next_idx[-1] - limited_idx_num * num
                        binary_sequence = bin(0)[2:].zfill(first_idx_length) \
                                          + 'xxxx' + bin(num+add_next_idx[0])[2:].zfill(next_idx_length - 4) \
                                          + 'xxxx' + bin(final_idx)[2:].zfill(next_idx_length - 4) \
                                          + binary_string[c * payload_length:(c + 1) * payload_length]
                        modified_binary_lst.append(binary_sequence)

        new_version_dnas = new_version_dnas + modified_binary_lst

    indexs_binary, data_set, indexes, = index_operator.divide_all_version(new_version_dnas, need_log,
                                                                          first_idx_length, next_idx_length)

    _, sort_binary_lst = index_operator.sort_order_version(indexes, data_set, indexs_binary, need_log)
    dna_sequences = method.encode(sort_binary_lst, need_log)
    data_handle.write_dna_file(output_path, dna_sequences, need_log)


def decode(method, new_version_dnas_path, output_path, model_path=None, need_index=True,
           first_idx_length=20, next_idx_length=14, need_log=False):
    """
    Use the selected method, convert DNA sequences to the original file.

    :param method: transcoding method.
                   If you have model file, you can use this function without method.
    :type: object

    :param new_version_dnas_path: the path of DNA sequence set you need to convert.
    :type: string

    :param output_path: the path of binary sequences file.
    :type: string

    :param model_path: the path of model file if you want to save.
    :type: string

    :param verify: error correction method.
    :type: object

    :param need_index: whether the DNA sequences contain binary sequence indexes.
    :type: bool

    :param first_idx_length: the length of the prime index.
    :type: int

    :param next_idx_length: the total length of the minor index.
                            Note that the length of the minor index is the length of the marking plus the valid length.
    :type: int

    :param need_log: show the log.
    :type: bool
    """
    _, dna_sequences = data_handle.read_dna_file(new_version_dnas_path, need_log)
    output_binary_lst, size = method.decode(dna_sequences, need_log)

    if need_index:
        indexs_binary, data_set, indexes = index_operator.divide_all_version(output_binary_lst, need_log,
                                                                             first_idx_length, next_idx_length)
        output_binary_lst, _ = index_operator.sort_order_version(indexes, data_set, indexs_binary, need_log)

    data_handle.write_all_from_binary(output_path, output_binary_lst, size, need_log)
