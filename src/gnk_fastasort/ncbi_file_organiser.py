import sys
import logging

from gnk_fastasort.fetch_reports import fetch_sequence_reports
from gnk_fastasort.generics import compute_new_order, ordered_list, compute_group_lengths, calc_total_length, generate_final_dict, flatten_unlocs

logger = logging.getLogger('gnk_fastasort_logger')

def convert_names(order_dict, index):
    """
    Convert name from SUPER_* style to Genbank accession to reorder the fasta
    """
    new_final_dict = {}
    for i in order_dict:
        new_final_dict[index[i]["genbank_accession"]] = { 'original': i,'parent': order_dict[i]["parent"], 'length': order_dict[i]["length"]}

    return new_final_dict

def shrink_the_report(accession, report):
    report_condensed = {}
    for chromosome in report:
        if chromosome["assembly_accession"] == accession:
            report_condensed[chromosome["sequence_name"]] = {
                'genbank_accession': chromosome['genbank_accession'],
                'length': chromosome['length'],
                'molecule_type': chromosome['assigned_molecule_location_type']
            }
        else:
            sys.exit(f"[shrink_the_report] Provided accession {accession} doesn't match accession in report {chromosome['assembly_accession']}")
    return report_condensed

def main(args):
    report = fetch_sequence_reports(args.gca_accession)
    report_condensed: dict  = shrink_the_report(args.gca_accession, report)

    logger.info(f"[fastasort] Report condensed, record count: {len(report_condensed)}")

    new_order = compute_new_order(report_condensed)
    output_list: list[str] = ordered_list(new_order)

    group_lengths = compute_group_lengths(new_order, report_condensed)
    total_genome_length: int = calc_total_length(group_lengths)

    logger.info(f"[fastasort] Total Genome Length: {total_genome_length/1e6:.2f} Mb")

    # Sort the chrom_groups by length
    groups_sorted = sorted(group_lengths.items(), key=lambda x: (-x[1], x[0]))
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
            for name, data in report_condensed.items()
            if name not in ordered_list(new_order)
    ]

    # Sort the minor scaffolds by length
    minor_scaffolds_sorted: list[tuple[str, int]] = sorted(minor_scaffolds, key=lambda x: (-x[1], x[0]))

    # Return just the names of the minor scaffolds
    minor_scaffolds_names = [name for name, _ in minor_scaffolds_sorted]
    logger.info(f"[fastasort] Sorted MINOR Scaffold count: {len(minor_scaffolds_names)}")
    logger.info(f"[fastasort] Total Scaffolds (recalculated): {len(minor_scaffolds_names) + len(groups_sorted)}")

    # Merge the sorted lists of names
    sorted_all_groups: list[tuple[str, int]] = output_list + minor_scaffolds_names

    # flatten a dict of dict+list (parent molecule = { unlocs = [item1...]})
    # into a dict of unloc = {parent molecule}
    flattened_unlocs = flatten_unlocs(new_order_sorted)

    # Format the dictionary into a informative format
    final_dict = generate_final_dict(sorted_all_groups, flattened_unlocs, report_condensed, minor_scaffolds_names)

    logger.info(f"[fastasort] Final scaffold count: {len(final_dict)}")

    renamed_final_dict = convert_names(final_dict, report_condensed)

    with open(f"{args.gca_accession}_reordered.tsv", "w") as output:
        for x, y in renamed_final_dict.items():
            if args.style == "names":
                output.write(f"{x}\n")
            else:
                output.write(f"{x}\t{y['original']}\t{y['parent']}\t{y['length']}\n")
