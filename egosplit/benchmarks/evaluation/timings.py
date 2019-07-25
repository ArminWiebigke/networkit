from cover_benchmark import CoverBenchmark
from evaluation.output import create_line


def write_timings(benchmarks, result_dir, append):
	open_mode = 'a' if append else 'w'
	out_file = open(result_dir + 'timings.result', open_mode)
	if not append:
		print_header(out_file)
	for bm in benchmarks:
		timings = bm.get_timings()
		for name, value in timings.items():
			name_str = str(name, 'utf-8').replace(" ", "_")
			out_file.write(create_line(*bm.output_line(), name_str, value))


# Print output file headers
def print_header(out_file):
	out_file.write(create_line(*CoverBenchmark.output_header(), 'Timer Name', 'Running Time'))
