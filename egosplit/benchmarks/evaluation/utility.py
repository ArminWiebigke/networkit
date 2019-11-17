def create_line(*args, sep=','):
	line = ''
	for arg in args:
		line += str(arg) + sep
	line = line[:-len(sep)]
	line += '\n'
	return line
