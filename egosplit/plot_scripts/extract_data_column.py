def create_new_column(data, name, create_paras):
	# Create new data columns for x or hue if necessary
	# Example: Create the x column by taking the value between '_f-' and the next '*'
	#   from the 'algo' column
	# x = {'name': 'factor', 'create_from': 'algo', 'str_start': '_f-', 'str_end': '*',
	#      'type': 'float'}
	for row in data.itertuples(index=True, name='Pandas'):
		index = row.Index
		from_string = data.loc[index, create_paras['create_from']]

		str_start = create_paras['str_start']
		start_idx = from_string.find(str_start)
		if start_idx == -1:
			continue

		start_idx += len(str_start)
		str_end = create_paras['str_end']
		end_idx = from_string.find(str_end, start_idx)
		if end_idx == -1 or end_idx == start_idx:
			end_idx = len(from_string)

		value = from_string[start_idx:end_idx]
		value = create_paras['type'](value)
		data.loc[index, name] = value
