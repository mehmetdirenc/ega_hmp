import os

def create_taxonomy_dict(tax_file):
    taxonomy_dict = {}
    with open(tax_file, 'r') as fh:
        for line in fh:
            if line.startswith('Feature'):
                continue
            id = line.split('\t')[0]
            tax = line.split('\t')[1].split(";")[-1].split("__")[1].strip()
            taxonomy_dict[id] = tax
    return taxonomy_dict

def change_newick_ids(tree_path, new_tree_path, tax_file):
    taxonomy_dict = create_taxonomy_dict(tax_file)
    with open(tree_path, 'r') as fh, open(new_tree_path, 'w') as fw:
        for line in fh:
            for tax_id in taxonomy_dict:
                line = line.replace(tax_id, taxonomy_dict[tax_id])
            fw.write(line)




if __name__ == '__main__':
    tax_file = "/mnt/lustre/home/mager/magmu818/datasets/public_databases/gtdb/release220/taxonomy/bac120_taxonomy_r220.tsv"
    tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode7/gtdbtk/tree/classify/itol_barcode7.tree"
    new_tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode7/gtdbtk/tree/classify/itol_barcode7_replaced.tree"
    change_newick_ids(tree_path, new_tree_path, tax_file)
    tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode9/gtdbtk/tree/classify/itol_barcode9.tree"
    new_tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode9/gtdbtk/tree/classify/itol_barcode9_replaced.tree"
    change_newick_ids(tree_path, new_tree_path, tax_file)
    tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode8/gtdbtk/tree/classify/itol_gtdbtk.bac120.classify.tree.2.tree"
    new_tree_path = "/mnt/lustre/home/mager/magmu818/results/mahana/WGS/hybrid/barcode8/gtdbtk/tree/classify/itol_barcode8_replaced.tree"
    change_newick_ids(tree_path, new_tree_path, tax_file)