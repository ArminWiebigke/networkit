class A:
	msg = "A"

	def print(self):
		print(self.msg)

	def set(self):
		msg = "X"


a = A()
a2 = A()
a.msg = "B"
a.print()
a2.print()