import sys
from Bio import Phylo
import os
import csv

EUROPE_COUNTRIES = set([
    "UK","United Kingdom","Norway","Netherlands","Switzerland",
    "Germany","France","Italy","Spain","Denmark","Sweden",
    "Finland","Belgium","Austria","Ireland","Portugal",
    "Poland","Czech Republic","Hungary","Greece"
])

def get_terminals_under_clade(clade):
    return list(clade.find_clades(terminal=True))

def load_metadata(metadata_file):
    metadata = {}
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            strain = row['Strain']
            country = row['Country']
            migration = row['Migration']
            metadata[strain] = (country, migration)
    return metadata

def select_nodes(tree, node_info, threshold):
    all_terminals = set(tree.get_terminals())
    selected_nodes = []
    covered_terminals = set()
    eligible_nodes = [(node, info)
                      for node, info in node_info.items()
                      if info[0] < threshold]
    sorted_nodes = sorted(
        eligible_nodes,
        key=lambda x: -x[1][2]
    )

    for node, _ in sorted_nodes:
        if covered_terminals == all_terminals:
            break
        node_terminals = set(
            get_terminals_under_clade(node)
        )
        if not node_terminals.issubset(covered_terminals):
            selected_nodes.append(node)
            covered_terminals.update(node_terminals)

    for terminal in all_terminals - covered_terminals:
        selected_nodes.append(terminal)
    return selected_nodes

def classify_cluster(terminals, metadata):

    if len(terminals) == 1:

        strain = terminals[0]

        if strain in metadata:

            country, migration = metadata[strain]

            return {
                strain: (migration, migration)
            }

        else:

            return {
                strain: ("Unknown","Unknown")
            }

    countries = []

    migrations = []

    for strain in terminals:

        if strain in metadata:

            country, migration = metadata[strain]

            countries.append(country)

            migrations.append(migration.lower())

    has_non_europe = any(
        c not in EUROPE_COUNTRIES
        for c in countries
    )

    has_immigrant = any(
        m == "immigrant"
        for m in migrations
    )


    updated_status = {}


    for strain in terminals:

        if strain not in metadata:

            updated_status[strain] = ("Unknown","Unknown")

            continue


        country, migration = metadata[strain]

        migration_lower = migration.lower()

        if migration_lower == "unknown":

            if has_non_europe or has_immigrant:

                updated = "Potential immigrant"

            else:

                updated = "Potential local"

        else:

            updated = migration


        updated_status[strain] = (
            migration,
            updated
        )

    return updated_status

def calculate_clusters(input_file,
                       threshold,
                       metadata_file):

    tree = Phylo.read(input_file,
                      "newick")

    metadata = load_metadata(
        metadata_file
    )

    node_info = {}

    for clade in tree.find_clades():

        max_distance = 0

        terminal_names = []

        terminals = get_terminals_under_clade(clade)

        for terminal in terminals:

            distance = tree.distance(
                clade,
                terminal
            )

            max_distance = max(
                max_distance,
                distance
            )

            terminal_names.append(
                terminal.name
            )


        node_info[clade] = (

            max_distance,

            terminal_names,

            len(terminal_names)

        )


    selected_nodes = select_nodes(

        tree,

        node_info,

        threshold

    )


    output_file = os.path.splitext(
        input_file
    )[0] + "_cluster_migration.csv"

    with open(output_file,
              'w',
              newline='') as f:

        writer = csv.writer(f)

        writer.writerow([

            "ClusterID",

            "Strain",

            "Country",

            "Migration",

            "MigrationUpdated"

        ])

        cluster_id = 1

        for node in selected_nodes:


            if node.is_terminal():

                terminals = [node.name]

            else:

                terminals = [

                    t.name

                    for t in

                    get_terminals_under_clade(node)

                ]

            updated_status = classify_cluster(

                terminals,

                metadata
            )

            for strain in terminals:

                if strain in metadata:

                    country = metadata[strain][0]

                else:

                    country = "Unknown"

                migration, updated = updated_status[strain]

                writer.writerow([

                    cluster_id,

                    strain,

                    country,

                    migration,

                    updated

                ])

            cluster_id += 1

if __name__ == "__main__":

    if len(sys.argv) != 4:

        print("Usage:")

        print("python cluster_migration.py tree.nwk threshold metadata.csv")

        sys.exit(1)

    tree_file = sys.argv[1]

    threshold = float(sys.argv[2])

    metadata_file = sys.argv[3]

    calculate_clusters(

        tree_file,

        threshold,

        metadata_file

    )
