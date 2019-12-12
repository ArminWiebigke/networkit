from multiprocessing import Process
from time import sleep

from egosplit.benchmarks.data_structures.timeout import TimeoutException


def f(time):
	sleep(time)


def run_with_limited_time(func, args, kwargs, time):
	"""Runs a function with time limit

	:param func: The function to run
	:param args: The functions args, given as tuple
	:param kwargs: The functions keywords, given as dict
	:param time: The time limit in seconds
	:return: True if the function ended successfully. False if it was terminated.
	"""
	p = Process(target=func, args=args, kwargs=kwargs)
	p.start()
	p.join(time)
	if p.is_alive():
		print("Timed out!")
		p.terminate()
		return False
		# raise TimeoutException("Not finished")
	return True
