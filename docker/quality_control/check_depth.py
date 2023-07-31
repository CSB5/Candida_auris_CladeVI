#!/usr/bin/env python
# ========================================================
# Check sequencing depth and perform subsampling if needed
# ========================================================

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--genome_size", dest="genome_size", type=int, help="Genome size", required=True)
parser.add_argument("--total_bases", dest="total_bases", type=int, help="Total number of bases", required=True)
parser.add_argument("--min_target_depth", dest="min_target_depth", type=int, help="Minimum target depth", required=True)
parser.add_argument("--max_target_depth", dest="max_target_depth", type=int, help="Maximum target depth", required=True)

args = parser.parse_args()

current_depth = args.total_bases / args.genome_size
output = f"Current depth: {current_depth}\\n"

if current_depth < args.min_target_depth:
    # TODO create a warning
    output += f"-1"
elif current_depth < args.max_target_depth:
    # TODO record that the depth is already in the range
     output += f"1"
else:
    # TODO record that the depth exceed the maximum target depth and need to be subsampled
    proportion = args.max_target_depth / current_depth
    output += f"{proportion}"

# print(repr(output))
print(output)