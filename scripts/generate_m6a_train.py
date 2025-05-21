import sys
import csv
from tqdm import tqdm
import argparse
import os
import pandas as pd
import numpy as np
import random
import re

def parse_arguments():
	parser = argparse.ArgumentParser(description="Generate m6A sites training data.")
	parser.add_argument("--input_file", required=True, type=str, help="Input file after clean_event.py process")
	parser.add_argument("--output_file", required=True, type=str, help="Output dir path")
	parser.add_argument("--ref_pos_file", required=True, type=str, help="Reference position file")
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

# 读取结合位点文件，生成一个包含匹配项的集合
def read_binding_sites(binding_file):
	binding_sites = set()
	with open(binding_file, 'r') as file:
		reader = csv.reader(file, delimiter='\t')
		next(reader)  # 跳过标题行
		for row in reader:
			transcript_id = row[0]
			transcript_pos = int(row[1])
			binding_sites.add((transcript_id, transcript_pos))
	return binding_sites

# 处理电流数据文件，进行匹配并输出结果
def process_current_data(flt_modules, binding_sites):
	lines = flt_modules
	matched_modules = []
	unmatched_modules = []
	module = []
	
	# 处理标题行
	header = lines[0].strip()
	matched_modules.append(header)
	unmatched_modules.append(header)

	# 创建进度条
	total_lines = len(lines)
	for line in tqdm(lines[1:], desc="Processing modules", total=total_lines-1):  # 跳过标题行
		line = line.strip()
		if line == "----":  # 模块之间的分隔符
			if module:
				# 获取当前模块的第5个条目的 contig 和 position
				contig, position = module[4].split('\t')[:2]  # 修改为第5个条目
				contig = re.match(r'^ENST\d+\.\d+', contig).group(0)
				position = int(position)
				
				# 判断是否匹配
				if (contig, position) in binding_sites:
					matched_modules += module + ["----"]
				else:
					unmatched_modules += module + ["----"]
			module = []
		else:
			module.append(line)
	
	# 处理最后一个模块
	if module:
		contig, position = module[4].split('\t')[:2]  # 修改为第5个条目
		position = int(position)
		
		if (contig, position) in binding_sites:
			matched_modules += module + ["----"]
		else:
			unmatched_modules += module + ["----"]
	
	return matched_modules, unmatched_modules

def merge_module(flt_modules):

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
	if numeric_df.shape[0] != 0:
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
	else:
		print("Error: The merged module is empty.")
		sys.exit()
	return merge_df

def extract_random_entries(df, num_neg):    
	if num_neg > len(df):
		print(f"Requested number of entries ({num_neg}) exceeds total entries ({len(df)}). Extracting all entries.")
		num_neg = len(df)
	sampled_df = df.sample(n=num_neg, random_state=random.randint(0, 10000))
	
	return sampled_df

def merge_and_shuffle(df1, df2, output_file):
	
	df_merged = pd.concat([df1, df2], ignore_index=True)
	
	# 打乱行的顺序
	df_merged = df_merged.sample(frac=1, random_state=42).reset_index(drop=True)
	
	# 保存合并后的文件，并保留一个标题行
	df_merged.to_csv(output_file, sep="\t", index=False)
	print(f"Files merged and shuffled. Output saved as {output_file}.")

def main():
	args = parse_arguments()
	input_file = args.input_file
	output_file = args.output_file
	ref_pos_file = args.ref_pos_file

	if not os.path.exists(input_file):
		print(f"Error: The file {input_file} does not exist.")
		return

	module_list = extract_modules(input_file)
	print(f"Module extraction complete.\n")

	flt_modules = filter_modules(module_list)
	# 读取结合位点文件
	binding_sites = read_binding_sites(ref_pos_file)
	# 处理电流数据文件并生成输出
	matched_modules, unmatched_modules = process_current_data(flt_modules, binding_sites)
	# split pos and neg sample
	match_module_df = merge_module(matched_modules)
	unmatch_module_df = merge_module(unmatched_modules)
	
	# random sample neg data set
	# sampled_unmatch_df = extract_random_entries(unmatch_module_df, num_neg)

	# remove col and add label

	match_module_df['label'] = "1"
	unmatch_module_df['label'] = "0"

	# merge pos and neg sample
	merge_and_shuffle(match_module_df, unmatch_module_df, output_file)

if __name__ == "__main__":
	main()