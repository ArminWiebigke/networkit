from networkit.stopwatch import Timer
from egosplit.benchmarks.run import start_benchmarks

t = Timer()
start_benchmarks()
print("\nTotal time for benchmarks: {}s".format(str(t.stop())[:6]))