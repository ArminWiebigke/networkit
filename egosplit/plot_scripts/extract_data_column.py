def create_new_columns(data, hue, x):
	# Create new data columns for x or hue if necessary
	# Example: Create the x column by taking the value between '_f-' and the next '*'
	#   from the 'algo' column
	# x = {'name': 'factor', 'create_from': 'algo', 'str_start': '_f-', 'str_end': '*',
	#      'type': 'float'}
	new_columns = []
	if type(x) != str and x:
		new_columns.append(x)
	if type(hue) != str and hue:
		new_columns.append(hue)
	for new_column in new_columns:
		remove_rows = []
		for row in data.itertuples(index=True, name='Pandas'):
			index = row.Index
			from_string = data.loc[index, new_column['create_from']]

			str_start = new_column['str_start']
			start_idx = from_string.find(str_start)
			if start_idx == -1:
				raise RuntimeError("Could not create column {} from string \"{}\"".format(
					new_column['name'], from_string
				))
			start_idx += len(str_start)
			str_end = new_column['str_end']
			end_idx = from_string.find(str_end, start_idx)
			if end_idx == -1 or end_idx == start_idx:
				end_idx = len(from_string)

			value = from_string[start_idx:end_idx]
			value = new_column['type'](value)
			data.loc[index, new_column['name']] = value

		for row in remove_rows:
			data.drop(row, inplace=True)
