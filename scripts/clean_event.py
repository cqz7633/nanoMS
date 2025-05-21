import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from time import perf_counter

def parse_arguments():
	parser = argparse.ArgumentParser(description='Merge entries at the same location and perform data cleaning.')
	parser.add_argument('--input', type=str, required=True, help='Input data file path.')
	parser.add_argument('--output', type=str, required=True, help='Output result file path.')
	parser.add_argument('--processes', type=int, default=12, help='Number of processes used.')
	return parser.parse_args()

def filter_and_average(series, std_multiplier=3):
	"""优化后的过滤和平均函数"""
	if len(series) == 1:
		return series.iloc[0]
	
	median = np.median(series)
	std = np.std(series)
	
	if np.isnan(std) or std == 0:
		return median
	
	mask = np.abs(series - median) <= std_multiplier * std
	filtered = series[mask]
	
	return filtered.mean() if len(filtered) > 0 else median

def process_group(group):
	"""优化后的分组处理函数"""
	agg_data = {
		'reference_kmer': group['reference_kmer'].mode()[0] if not group['reference_kmer'].mode().empty else group['reference_kmer'].iloc[0],
		'event_level_mean': filter_and_average(group['event_level_mean']),
		'event_stdv': filter_and_average(group['event_stdv']),
		'event_length': filter_and_average(group['event_length']),
		'model_stdv': group['model_stdv'].iloc[0],
		'standardized_level': filter_and_average(group['standardized_level']),
		'event_level_median': group['event_level_mean'].median(),
		'event_stdv_median': group['event_stdv'].median(),
		'event_length_median': group['event_length'].median()
	}
	return agg_data

def process_group_wrapper(args):
	"""多进程包装函数"""
	name, group = args
	result = process_group(group)
	result['contig'], result['position'] = name
	return result

def main():
	args = parse_arguments()
	input_file = args.input
	output_file = args.output
	processes = args.processes

	# 读取数据
	print("Loading data...")
	try:
		df = pd.read_csv(input_file, sep='\t')
	except Exception as e:
		print(f"Error reading input file: {e}")
		return

	# remove col
	rm_col = [col for col in ["strand", "model_kmer", "model_mean"] if col in df.columns]
	df = df.drop(columns=rm_col)

	# 检查必要列
	required_columns = [
		'contig', 'position', 'reference_kmer',
		'event_level_mean', 'event_stdv', 'event_length', 'model_stdv', 'standardized_level'
	]
	missing_columns = [col for col in required_columns if col not in df.columns]
	if missing_columns:
		print(f"Lack of necessary columns: {missing_columns}")
		return

	# 分组处理
	print("Start grouping processing ...")
	grouped = df.groupby(['contig', 'position'])
	groups = list(grouped)
	
	# 多进程处理
	with Pool(processes=processes) as pool:
		results = list(tqdm(pool.imap(process_group_wrapper, groups), 
						  total=len(groups),
						  desc="Process grouping"))

	# 合并结果并调整列顺序
	print("Merge the results and save ...")
	result_df = pd.DataFrame(results)
	
	# 确保contig和position列在最前面
	columns_order = ['contig', 'position'] + [col for col in result_df.columns if col not in ['contig', 'position']]
	result_df = result_df[columns_order]
	
	### clean1
	result_df["diff_stdv"] = result_df["event_stdv"] - result_df["model_stdv"]
	result_df["diff_stdv"] = result_df["diff_stdv"] / result_df["model_stdv"]
	
	### clean2
	if 'model_stdv' in result_df.columns:
		result_df = result_df.drop(columns=['model_stdv'])
	# result_df['contig'] = result_df['contig'].apply(lambda x: x.split('|')[0])
	
	result_df.to_csv(output_file, sep='\t', index=False)
	print(f"结果已保存到 {output_file}")

if __name__ == "__main__":
	main()

# import pandas as pd
# import numpy as np
# import argparse
# from tqdm import tqdm
# from multiprocessing import Pool, cpu_count
# import os
# from time import perf_counter
# import psutil
# def parse_arguments():
#     parser = argparse.ArgumentParser(description='合并同一位置的条目，并进行数据清洗。')
#     parser.add_argument('--input', type=str, required=True, help='已排序的输入数据文件路径（TSV 格式）')
#     parser.add_argument('--output', type=str, required=True, help='输出结果文件路径（TSV 格式）')
#     parser.add_argument('--processes', type=int, default=20, help='使用的进程数')
#     parser.add_argument('--chunk-size', type=int, default=10000000, help='每个数据块的行数')
#     return parser.parse_args()

# def filter_and_average(series, std_multiplier=3):
#     """优化后的过滤和平均函数"""
#     if len(series) == 1:
#         return series.iloc[0]
	
#     median = np.median(series)
#     std = np.std(series)
	
#     if np.isnan(std) or std == 0:
#         return median
	
#     mask = np.abs(series - median) <= std_multiplier * std
#     filtered = series[mask]
	
#     return filtered.mean() if len(filtered) > 0 else median

# def process_group(group):
#     """优化后的分组处理函数"""
#     agg_data = {
#         'reference_kmer': group['reference_kmer'].mode()[0] if not group['reference_kmer'].mode().empty else group['reference_kmer'].iloc[0],
#         'strand': group['strand'].iloc[0],
#         'event_level_mean': filter_and_average(group['event_level_mean']),
#         'event_stdv': filter_and_average(group['event_stdv']),
#         'event_length': filter_and_average(group['event_length']),
#         'model_kmer': group['model_kmer'].iloc[0],
#         'model_mean': group['model_mean'].iloc[0],
#         'model_stdv': group['model_stdv'].iloc[0],
#         'standardized_level': filter_and_average(group['standardized_level']),
#         'event_level_median': group['event_level_mean'].median(),
#         'event_stdv_median': group['event_stdv'].median(),
#         'event_length_median': group['event_length'].median()
#     }
#     return agg_data

# def process_group_wrapper(args):
#     """多进程包装函数"""
#     name, group = args
#     result = process_group(group)
#     result['contig'], result['position'] = name
#     return result

# def get_memory_usage():
#     process = psutil.Process(os.getpid())
#     return process.memory_info().rss / 1024**2  # 返回MB

# def main():
#     args = parse_arguments()
	
#     # 检查必要列
#     required_columns = ['contig', 'position', 'reference_kmer', 'strand',
#                       'event_level_mean', 'event_stdv', 'event_length',
#                       'model_kmer', 'model_mean', 'model_stdv', 'standardized_level']
	
#     # 初始化缓冲区
#     buffer = pd.DataFrame()
#     first_chunk = True
	
#     # 读取数据块
#     for chunk in pd.read_csv(args.input, sep='\t', chunksize=args.chunk_size):
#         # 合并缓冲区数据
#         if not buffer.empty:
#             chunk = pd.concat([buffer, chunk], ignore_index=True)
		
#         # 处理空数据块
#         if chunk.empty:
#             continue
		
#         # 查找最后一个完整组的位置
#         last_contig = chunk['contig'].iloc[-1]
#         last_position = chunk['position'].iloc[-1]
		
#         # 获取最后组的所有行
#         mask = (chunk['contig'] == last_contig) & (chunk['position'] == last_position)
#         last_group_indices = chunk[mask].index
		
#         if not last_group_indices.empty:
#             last_group_end = last_group_indices[-1] + 1
#             process_part = chunk.iloc[:last_group_end]
#             buffer = chunk.iloc[last_group_end:]
#         else:
#             process_part = chunk
#             buffer = pd.DataFrame()
		
#         # 处理当前数据块
#         if not process_part.empty:
#             # 检查列是否完整
#             missing_columns = [col for col in required_columns if col not in process_part.columns]
#             if missing_columns:
#                 print(f"缺少必要的列: {missing_columns}")
#                 return
			
#             # 分组处理
#             grouped = process_part.groupby(['contig', 'position'])
#             groups = list(grouped)
			
#             # 多进程处理
#             with Pool(processes=args.processes) as pool:
#                 results = list(tqdm(pool.imap(process_group_wrapper, groups), 
#                                   total=len(groups),
#                                   desc="处理分组"))
			
#             # 合并结果并调整列顺序
#             result_df = pd.DataFrame(results)
#             columns_order = ['contig', 'position'] + [col for col in result_df.columns if col not in ['contig', 'position']]
#             result_df = result_df[columns_order]
#             print(f"处理内存: {get_memory_usage():.2f} MB")
#             # 写入输出文件
#             if first_chunk:
#                 result_df.to_csv(args.output, sep='\t', index=False)
#                 first_chunk = False
#             else:
#                 result_df.to_csv(args.output, sep='\t', mode='a', header=False, index=False)
	
#     # 处理缓冲区剩余数据
#     if not buffer.empty:
#         # 检查列是否完整
#         missing_columns = [col for col in required_columns if col not in buffer.columns]
#         if missing_columns:
#             print(f"缺少必要的列: {missing_columns}")
#             return
		
#         # 分组处理
#         grouped = buffer.groupby(['contig', 'position'])
#         groups = list(grouped)
		
#         # 多进程处理
#         with Pool(processes=args.processes) as pool:
#             results = list(tqdm(pool.imap(process_group_wrapper, groups), 
#                               total=len(groups),
#                               desc="处理剩余分组"))
		
#         # 合并结果并调整列顺序
#         result_df = pd.DataFrame(results)
#         columns_order = ['contig', 'position'] + [col for col in result_df.columns if col not in ['contig', 'position']]
#         result_df = result_df[columns_order]
		
#         # 写入输出文件
#         result_df.to_csv(args.output, sep='\t', mode='a', header=False, index=False)

# if __name__ == "__main__":
#     start = perf_counter()
#     main()
#     print(f"总耗时：{perf_counter() - start:.6f} 秒")