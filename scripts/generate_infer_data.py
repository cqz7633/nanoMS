import sys
import csv
from tqdm import tqdm
import argparse
import os
import pandas as pd
import numpy as np

def parse_arguments():
	parser = argparse.ArgumentParser(description="Filter modules based on DRACH rule.")
	parser.add_argument("--input_file", type=str, required=True, help="Input file after clean_event.py process")
	parser.add_argument("--output_dir", type=str, required=True, help="Output dir path")
	parser.add_argument("--file_prefix", type=str, required=True, help="Output file prefix")
	return parser.parse_args()

def extract_modules(input_file):
	module_list = []
	with open(input_file, 'r') as infile:
		reader = csv.DictReader(infile, delimiter='\t')
		#writer = csv.writer(outfile, delimiter='\t')

		module_list.append("\t".join(reader.fieldnames))
		rows = list(reader)
		total_rows = len(rows)
		
		progress = tqdm(range(3, total_rows - 3), desc="Processing", unit="module")
		for i in progress:
			module = rows[i-3:i+4]
			contigs = {row['contig'] for row in module}
			if len(contigs) != 1:
				continue

			
			positions = [int(row['position']) for row in module]
			expected_positions = list(range(positions[0], positions[0] + 7))
			if positions != expected_positions:
				continue
			for row in module:
				module_list.append("\t".join([row[field] for field in reader.fieldnames]))
			module_list.append("----")

		progress.close()
	return module_list

valid_first_bases = {
	'D': ['A', 'G', 'T'],  # D -> A, G, T (T is used for U)
	'R': ['A', 'G'],       # R -> A, G
	'A': ['A'],            # A -> A
	'C': ['C'],            # C -> C
	'H': ['A', 'C', 'T']   # H -> A, C, T (T is used for U)
}

# 判断一个碱基是否符合DRACH规则
def is_valid_drach_sequence(kmer_sequence):

	if len(kmer_sequence) != 5:
		return False  # 确保长度为5

	# 第一个碱基：符合D规则
	if kmer_sequence[0] not in valid_first_bases['D']:
		return False
	# 第二个碱基：符合R规则
	if kmer_sequence[1] not in valid_first_bases['R']:
		return False
	# 第三个碱基：必须是A
	if kmer_sequence[2] != 'A':
		return False
	# 第四个碱基：必须是C
	if kmer_sequence[3] != 'C':
		return False
	# 第五个碱基：符合H规则
	if kmer_sequence[4] not in valid_first_bases['H']:
		return False
	
	return True

def filter_modules(module_list):
	flt_modules = []
	
	modules = []
	current_module = []
	
	for line in module_list:

		if line.startswith("----"):  
			if current_module:
				modules.append(current_module)
			current_module = []
		elif line:  # 非空行
			current_module.append(line.split())


	if current_module:
		modules.append(current_module)

	# f.write("\t".join(module_list[0].split()) + "\n")
	flt_modules.append("\t".join(module_list[0].split()))
	for module in tqdm(modules, desc="Processing modules", unit="module"):
		kmer_sequence = module[1][2]
		#kmer_sequence = ''.join([row[2][0] for row in module[1:6]])  
		if is_valid_drach_sequence(kmer_sequence):
			for row in module:
				flt_modules.append("\t".join(row))
				# f.write("\t".join(row) + "\n")
			# f.write("----\n")
			flt_modules.append("----")
	return flt_modules

def merge_module(flt_modules, output_file):

	header = flt_modules[0].split('\t')

	modules = []
	current_module = []
	for line in flt_modules[1:]:
		# line = line.strip()
		if line == '----':
			if current_module:
				modules.append(current_module)
			current_module = []
		elif line:
			current_module.append(line.split('\t'))
	if current_module:
		modules.append(current_module)

	new_header = ["contig", "position"]
	data_cols = header[2:]
	offsets = range(-3, 4)
	for i in offsets:
		prefix = f"+{i}" if i > 0 else str(i)
		for col_name in data_cols:
			new_header.append(f"{prefix}_{col_name}")
	merge_list = []
	for module in tqdm(modules, desc="Processing Modules"):
		fourth_line = module[3]
		contig_4th = fourth_line[0]
		position_4th = fourth_line[1]
		merged_row = [contig_4th, position_4th]
		for idx, row in enumerate(module):
			data_values = row[2:]
			merged_row.extend(data_values)
		merge_list.append(merged_row)
		# out_file.write('\t'.join(merged_row) + '\n')
	merge_df = pd.DataFrame(merge_list, columns = new_header)
	columns_to_drop = [col for col in merge_df.columns if 'reference_kmer' in col]
	merge_df = merge_df.drop(columns=columns_to_drop)

	# filter nan inf
	numeric_df = merge_df[merge_df.columns[2:]].apply(pd.to_numeric, errors='coerce')
	mask = (
		numeric_df.isna().any(axis=1) |          # 检查NaN
		np.isinf(numeric_df).any(axis=1)        # 检查inf/-inf
	)
	merge_df = merge_df[~mask]
	removed_count = mask.sum()
	if removed_count > 0:
		print(f"Remove {removed_count} lines containing NaN or inf/-inf.")
	else:
		print("No samples containing NaN or nf/- nf were found.")
	merge_df["label"] = -1
	merge_df.to_csv(output_file, sep="\t", index=False)

def main():
	args = parse_arguments()
	input_file = args.input_file
	output_dir = args.output_dir
	output_prefix = args.file_prefix

	if not os.path.exists(input_file):
		print(f"Error: The file {input_file} does not exist.")
		return

	if not os.path.exists(output_dir):
		print(f"Error: Path '{output_dir}' does not exist！")
		sys.exit(1) 
	site_output = os.path.join(output_dir, output_prefix+".m6a.tsv") 
	stru_output = os.path.join(output_dir, output_prefix+".structure.tsv") 

	module_list = extract_modules(input_file)
	print(f"Module extraction complete.\n")

	merge_module(module_list, stru_output)
	print(f"Merge structure modules saved to {stru_output}.")

	flt_modules = filter_modules(module_list)
	print(f"Filtered DRACH sites.\n")

	merge_module(flt_modules, site_output)
	print(f"Merge m6a sites modules saved to {site_output}.")

if __name__ == "__main__":
	main()