def create_line(*args, sep=','):
	line = ''
	for arg in args:
		if isinstance(arg, str) and sep in arg:
			raise RuntimeError("Seperator '{}' used in value!".format(sep))
		line += str(arg) + sep
	line = line[:-len(sep)]
	line += '\n'
	return line
