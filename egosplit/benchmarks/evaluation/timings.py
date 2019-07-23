from evaluation.output import create_line


def write_timings(benchmarks, result_dir, append, file_suffix):
	open_mode = 'a' if append else 'w'
	out_file = open(result_dir + 'timings.result', open_mode)
	if not append:
		print_header(out_file)
	for benchmark in benchmarks:
		timings = benchmark.get_timings()
		for name, value in timings.items():
			name_str = str(name, 'utf-8').replace(" ", "_")
			out_file.write(create_line(benchmark.get_algo_name(), benchmark.get_graph_name(),
			                           benchmark.get_graph_id(), name_str, value))


# Print output file headers
def print_header(out_file):
	out_file.write(create_line('algo', 'graph', 'graph_id', 'timer_name', 'time'))
