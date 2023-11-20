from yyc import pipeline
from yyc import scheme
from yyc.utils import data_handle


# ## encode original version file
# pipeline.encode_original(method=scheme.YYC(support_bases="A", base_reference=[0, 1, 0, 1],
#                                            current_code_matrix=[[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]],
#                                            search_count=100, max_homopolymer=4, max_content=0.6),
#                          input_path="./examples/files/version1.txt",
#                          output_path="./examples/output/version1.dna",
#                          model_path="./examples/output/version1.pkl",
#                          need_index=True, need_log=True)
# output_version = 1
# output_file_path = "./examples/output/version" + str(output_version) + ".txt"

# # version2
# operation = ['insert'] + ['replace']*2
# last_version_dna = "version1.dna"
# output_version = 2

# # version3
# operation = ['insert']
# last_version_dna = "version2.dna"
# output_version = 3

# version4
operation = ['replace']*4 + ['insert']
last_version_dna = "version3.dna"
output_version = 4

root_path = "./examples/files/match_files/version" + str(output_version) + "/"
output_file_path = "./examples/output/version" + str(output_version) + ".txt"

## encode different version
pipeline.encode(method=scheme.YYC(support_bases="A", base_reference=[0, 1, 0, 1],
                                  current_code_matrix=[[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]],
                                  search_count=100, max_homopolymer=4, max_content=0.6),
                modifications=operation,
                modified_files=[root_path + "v" + str(output_version) + "_" + str(i+1) + "_" + operation[i] + ".txt"
                                if operation[i] != 'replace'
                                else [root_path + "v" + str(output_version) + "_" + str(i+1) + "_replace_delete.txt",
                                      root_path + "v" + str(output_version) + "_" + str(i+1) + "_replace_insert.txt"]
                                for i in range(len(operation))],
                match_files=[root_path + "match_v" + str(output_version) + "_" + str(i+1) + ".txt"
                             for i in range(len(operation))],
                last_version_dnas_path="./examples/output/" + last_version_dna,
                output_path="./examples/output/" + "version" + str(output_version) + ".dna")


## decode
pipeline.decode(method=scheme.YYC(support_bases="A", base_reference=[0, 1, 0, 1],
                                  current_code_matrix=[[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]],
                                  search_count=100, max_homopolymer=4, max_content=0.6),
                new_version_dnas_path="./examples/output/version" + str(output_version) + ".dna",
                output_path=output_file_path, need_index=True, need_log=True)
