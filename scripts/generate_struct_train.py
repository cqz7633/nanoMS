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
	parser = argparse.ArgumentParser(description="Filter modules based on DRACH rule.")
	parser.add_argument("--input_file", required=True, type=str, help="Input file after clean_event.py process")
	parser.add_argument("--output_file", required=True, type=str, help="Output file path")
	parser.add_argument("--shape_file", required=True, type=str, help="Reference icshape file")
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


# 读取结合位点文件，生成一个包含匹配项的集合
def read_icshape(icshape_file):
	struct_dict = {}
	with open(icshape_file, "r") as sf:
		for line in sf:
			line = line.strip()
			if not line:
				continue
			cols = line.split('\t')
			if len(cols) < 4:
				contig = cols[0]
				struct_dict[contig] = []
				continue
			contig = cols[0]
			structure_scores = cols[3:]
			if contig not in struct_dict:
				struct_dict[contig] = structure_scores
	return struct_dict

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
				print(contig)
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

	header = "\t".join(merge_df.columns)

	# 2. 用 "\t".join() 拼接每一行
	rows = ["\t".join(map(str, row)) for row in merge_df.values]

	# 3. 组合成列表，列名作为第一个元素
	res_list = [header] + rows
	return res_list

def extract_shape(merge_list, struct_dict):
	add_struct_list = []
	lines = merge_list
	if not lines:
		print("The current file is empty and has been exited.")
		return
	header = lines[0].split('\t')
	header.append("label")
	add_struct_list.append("\t".join(header))

	for line in tqdm(lines[1:], desc="Processing current file"):
		line = line.strip()
		if not line:
			continue
		cols = line.split('\t')
		if len(cols) < 2:
			cols.append("NULL")
			add_struct_list.append("\t".join(cols))
			continue
		contig = cols[0]
		contig = re.match(r'^ENST\d+\.\d+', contig).group(0)
		pos_str = cols[1]
		try:
			position = int(pos_str)
		except ValueError:
			cols.append("NULL")
			add_struct_list.append("\t".join(cols))
			continue
		label = "NULL"
		if contig in struct_dict:
			scores = struct_dict[contig]
			idx = position
			if 0 <= idx < len(scores):
				label = scores[idx]
		cols.append(label)
		add_struct_list.append("\t".join(cols))
	return add_struct_list

def is_valid_label(value):
	if value == "NULL":
		return False
	try:
		num = float(value)
		return num <= 0.2 or num >= 0.8
	except ValueError:
		return False

def filter_structure(add_struct_list):
	lines = add_struct_list
	if not lines:
		print("The input file is empty.")
		return

	header = lines[0]
	data_lines = lines[1:]
	flt_struct_list = []
	flt_struct_list.append(header)
	for line in tqdm(data_lines, desc="Filtering lines"):
		if not line.strip():
			continue
		cols = line.split('\t')
		label = cols[-1]  # 最后一列是 label
		if is_valid_label(label):
			flt_struct_list.append(line)
	return flt_struct_list

def convert_label(label_value):
	if label_value == "NULL":
		return label_value
	try:
		num = float(label_value)
		if num <= 0.2:
			return "0"
		elif num >= 0.8:
			return "1"
		else:
			return label_value
	except ValueError:
		return label_value

def convert_structure(flt_struct_list, output_file):
	lines = flt_struct_list
	if not lines:
		print("The input file is empty.")
		return

	header = lines[0]
	data_lines = lines[1:]
	conver_struct_list = []
	conver_struct_list.append(header)
	for line in tqdm(data_lines, desc="Converting label values"):
		if not line.strip():
			continue
		cols = line.split('\t')
		# 处理最后一列（label）
		cols[-1] = convert_label(cols[-1])
		conver_struct_list.append('\t'.join(cols))
	with open(output_file, 'w') as out:
		for line in tqdm(conver_struct_list, desc="Writing output"):
			out.write(line + '\n')

# def remove_columns(conver_struct_list, output_file):
# 	lines = conver_struct_list
# 	if not lines:
# 		print("The input file is empty.")
# 		return

# 	header = lines[0].split('\t')

# 	# 找到要移除的列的索引
# 	remove_cols = {'contig', 'position'}
# 	remove_indices = [i for i, col in enumerate(header) if col in remove_cols]

# 	# 构造新表头
# 	new_header = [col for i, col in enumerate(header) if i not in remove_indices]

# 	new_data_list = []
# 	new_data_list.append('\t'.join(new_header))

# 	for line in tqdm(lines[1:], desc="Processing lines"):
# 		line = line.strip()
# 		if not line:
# 			continue
# 		cols = line.split('\t')
# 		new_cols = [val for i, val in enumerate(cols) if i not in remove_indices]
# 		new_data_list.append('\t'.join(new_cols))

#     with open(output_file, 'w') as out:
# 		for line in tqdm(new_data_list, desc="Writing output"):
# 			out.write(line + '\n')

def balance_labels(new_data_list, output_file, seed=42):
	random.seed(seed)

	lines = new_data_list
	if not lines:
		print("The input file is empty.")
		return

	header = lines[0]
	data_lines = lines[1:]

	label_1 = []
	label_0 = []

	for line in tqdm(data_lines, desc="Parsing lines"):
		line = line.strip()
		if not line:
			continue
		cols = line.split('\t')
		label = cols[-1]
		if label == "1":
			label_1.append(line)
		elif label == "0":
			label_0.append(line)

	n_1 = len(label_1)
	print(f"Label=1 entries: {n_1}")

	if len(label_0) < n_1:
		print(f"⚠️ Warning: There are less than {n_1} entries with label=0, only {len(label_0)} entries, all of which will be used.")
		sampled_label_0 = label_0
	else:
		sampled_label_0 = random.sample(label_0, n_1)

	print(f"Randomly selected label=0 entries: {len(sampled_label_0)}")

	# 合并 + 打乱（如不需要可去除）
	final_data = label_1 + sampled_label_0
	random.shuffle(final_data)

	with open(output_file, 'w') as out:
		out.write(header + '\n')
		for line in tqdm(final_data, desc="Writing output"):
			out.write(line + '\n')

def main():

	args = parse_arguments()
	input_file = args.input_file
	output_file = args.output_file
	shape_file = args.shape_file


	if not os.path.exists(input_file):
		print(f"Error: The file {input_file} does not exist.")
		return

	module_list = extract_modules(input_file)
	print(f"Module extraction complete.\n")

	merge_list = merge_module(module_list)
	print(f"Merge structure modules.\n")

	struct_dict = read_icshape(shape_file)
	add_struct_list = extract_shape(merge_list, struct_dict)

	flt_struct_list = filter_structure(add_struct_list)
	
	convert_structure(flt_struct_list, output_file)
	
	# remove_columns(conver_struct_list, output_file)
	
	# balance_labels(new_data_list, output_file)

if __name__ == "__main__":
	main()