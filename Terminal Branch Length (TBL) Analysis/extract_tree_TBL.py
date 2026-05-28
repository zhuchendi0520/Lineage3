import sys
import csv
import re

def extract_tip_node_heights(tree_str):
    tip_node_heights = {}
    tip_pattern = re.compile(r'([\w.-]+):([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
    for match in tip_pattern.finditer(tree_str):
        tip_name, height = match.groups()
        if not tip_name.startswith("Node_"):  
            tip_node_heights[tip_name] = float(height)
    return tip_node_heights

if len(sys.argv) != 2:
    print("Usage: python script_name.py input_tree_file")
    sys.exit(1)

input_file = sys.argv[1]

with open(input_file, "r") as file:
    tree_str = file.read()

tip_node_heights = extract_tip_node_heights(tree_str)

output_file = f"{input_file.removesuffix('.tree')}_TBL.csv"

with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Tip Name", "Node Height"])
    for tip_name, height in tip_node_heights.items():
        writer.writerow([tip_name, height])
