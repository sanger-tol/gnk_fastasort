import sys
import logging

from gnk_fastasort.generics import (
    read_index, compute_new_order,
    ordered_list, compute_group_lengths,
    calc_total_length, generate_final_dict,
    flatten_unlocs
)

logger = logging.getLogger('gnk_fastasort_logger')

def main(args):
    if args.output:
        outname = f"{args.output}_reordered.tsv"
    else:
        outname = f"{str(args.index).split('.')[0]}_reordered.tsv"

    index: dict[str, dict[str, str]] = read_index(args.index)

    if len(index) == 0:
        logger.critical("[fastasort] INDEX IS EMPTY")
        sys.exit("[fastasort] INDEX IS EMPTY")

    logger.info(f"[fastasort] Index loaded, contains {len(index)} records")

    new_order: dict = compute_new_order(index)
    output_list: list[str] = ordered_list(new_order)

    group_lengths = compute_group_lengths(new_order, index)
    total_genome_length: int = calc_total_length(group_lengths)

    # Sort the chrom_groups by length
    groups_sorted = sorted(group_lengths.items(), key=lambda x: (-x[1], x[0]))

    logger.info(f"[fastasort] Sorted MAJOR Scaffold count: {len(output_list)}")
    logger.info(f"[fastasort] Total Genome Length: {total_genome_length/1e6:.2f} Mb")

    # Sort the chrom_groups by length
    groups_sorted: list[tuple[str, int]] = sorted(group_lengths.items(), key=lambda x: (-x[1], x[0]))
    logger.info(f"[fastasort] Sorted MAJOR groups: {groups_sorted}")

    # Ensure that the major groups are sorted as needed
    order_index = {name: idx for idx, (name, _) in enumerate(groups_sorted)}
    new_order_sorted = dict(
        sorted(new_order.items(), key=lambda kv: order_index.get(kv[0], float("inf")))
    )

    # non-SUPER named records
    # convert to a tuple of name, length
    minor_scaffolds: list[tuple[str, int]] = [
            (name, int(data["length"]))
            for name, data in index.items()
            if name not in ordered_list(new_order)
    ]

    # Sort the minor scaffolds by length
    minor_scaffolds_sorted: list[tuple[str, int]] = sorted(minor_scaffolds, key=lambda x: (-x[1], x[0]))

    # Return just the names of the minor scaffolds
    minor_scaffolds_names = [name for name, _ in minor_scaffolds_sorted]
    logger.info(f"[fastasort] Sorted MINOR Scaffold count: {len(minor_scaffolds_names)}")

    # Merge the sorted lists of names
    sorted_all_groups: list[str] = output_list + minor_scaffolds_names

    # flatten a dict of dict+list (parent molecule = { unlocs = [item1...]})
    # into a dict of unloc = {parent molecule}
    flattened_unlocs = flatten_unlocs(new_order_sorted)

    # Format the dictionary into a informative format
    final_dict = generate_final_dict(sorted_all_groups, flattened_unlocs, index, minor_scaffolds_names)

    logger.info(f"[fastasort] Final scaffold count: {len(final_dict)}")

    with open(outname, "w") as output:
        for i in final_dict:
            if args.style == "names":
                output.write(f"{i}\n")
            else:
                output.write(f"{i}\t-\t{final_dict[i]['parent']}\t{final_dict[i]['length']}\n")
