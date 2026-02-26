import re
import pathlib

def read_index(index_file: pathlib.Path) -> dict[str, dict[str, str]]:
    """
    Read the samtools faidx index file and return a dictionary
    of scaffold = length, offset, linebases, linewidth

    Standard: https://www.htslib.org/doc/faidx.html
    """
    with open(index_file, 'r') as f:
        index = {}
        for line in f:
            scaffold,  length, offset, linebases, linewidth = line.strip().split('\t')
            index[scaffold] = {
                "length": length,
                "offset": offset,
                "linebases": linebases,
                "linewidth": linewidth
            }
    return index

def compute_new_order(index: dict) -> dict[str, dict]:
    """
    Compute the new order of the scaffolds based on the
    naming of the scaffolds.
    """
    UNLOC_RE = re.compile(r'^(?P<prefix>[^\s>]+)_unloc_\d+$')
    scaff_groups: dict[str, dict] = {}

    for sequence_header in index.keys():
        unloc_match = UNLOC_RE.match(sequence_header)
        if unloc_match:
            prefix = unloc_match.group("prefix")
            scaff_groups.setdefault(prefix, {"main": None, "unlocs": []})
            scaff_groups[prefix]["unlocs"].append(sequence_header)
        elif sequence_header.startswith("SUPER_"):
            scaff_groups.setdefault(sequence_header, {"main": None, "unlocs": []})
            scaff_groups[sequence_header]["main"] = sequence_header
        else:
            pass

    return scaff_groups

def compute_group_lengths(scaff_group: dict, index: dict) -> dict:
    """
    Calculate the total length of each major group in dict
    """
    scaff_group_lengths: dict = {}
    for i in scaff_group:
        scaff_group_lengths[i] = int(index[i]["length"])
        for ii in scaff_group[i]["unlocs"]:
            scaff_group_lengths[i] += int(index[ii]["length"])

    return scaff_group_lengths

def calc_total_length(group_lengths: dict) -> int:
    """
    Sum the lenghts of all scaffolds
    """
    return sum(group_lengths.values())

def ordered_list(scaff_dictionary: dict) -> list[str] :
    """
    Return a list of scaffold names from a nested dict
    """
    ordered_list: list[str] = []
    for prefix, group in scaff_dictionary.items():
        ordered_list.append(prefix)
        ordered_list.extend(group["unlocs"])
    return ordered_list

def generate_final_dict(
    sorted_all_groups: dict,
    flattened_unlocs: dict,
    index: dict,
    minor_scaffolds_names: list[str]
) -> dict:
    """
    Generate a final dictionary with the parent and length of each scaffold
    """
    final_dict = {}

    for i in sorted_all_groups:
        if i in flattened_unlocs.keys():
            final_dict[i] = {
                "parent": flattened_unlocs[i]["parent"],
                "length": index[i]["length"]
            }
        elif i in index.keys() and i not in minor_scaffolds_names:
            final_dict[i] = {
                "parent": "MAJOR_UNIT",
                "length": index[i]["length"]
            }
        else:
            final_dict[i] = {
                "parent": "UNLOC_UNIT",
                "length": index[i]["length"]
            }

    return final_dict

def flatten_unlocs(data: dict) -> dict:
    """
    Convert:
      {"SUPER_14": {"main": "SUPER_14", "unlocs": ["SUPER_14_unloc_1"]}}
    into:
      {"SUPER_14_unloc_1": {"parent": "SUPER_14"}}
    """
    out = {}
    for parent, entry in data.items():
        for unloc in entry.get("unlocs", []):
            out[unloc] = {"parent": parent}
    return out
