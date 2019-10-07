from networkit.stopwatch import Timer


class ContextTimer:
	""" Time code like this:
	t = ContextTimer()
	with t:
		<Code>
	print(t.elapsed)
	"""

	def __init__(self):
		self.elapsed = 0.0
		self.timer = None

	def __enter__(self):
		self.timer = Timer()

	# On context exit, add the elapsed time since context enter
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.elapsed += self.timer.stop()
		del self.timer

	def reset(self):
		self.elapsed = 0.0
