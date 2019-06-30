/*
 * MemoizationTable.h
 *
 * Created: 2019-06-23
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_MEMOIZATIONTABLE_H
#define NETWORKIT_MEMOIZATIONTABLE_H

#include <vector>
#include <functional>

#include "../Globals.h"

namespace NetworKit {

template<typename T>
class MemoizationTable {
public:
	MemoizationTable(std::function<T(index)> valueFunc, T emptyVal, size_t initSize = 0);
	explicit MemoizationTable(T emptyVal, size_t initSize = 0);

	T getValue(index i);

	/**
	 * Try to set the function that should be memoized. This can only be done once.
	 * @param function
	 * @return true if function was set, false if the function was already set
	 */
	bool trySetValueFunc(std::function<T(index)> function);

private:
	std::function<T(index)> valueFunc;
	bool valueFuncSet;
	std::vector<T> data;
	T emptyVal;
};


template<typename T>
MemoizationTable<T>::MemoizationTable(std::function<T(index)> valueFunc, T emptyVal,
                                      size_t initSize)
		: valueFunc(valueFunc), valueFuncSet(true), emptyVal(emptyVal) {
	data.resize(initSize, emptyVal);
}

template<typename T>
MemoizationTable<T>::MemoizationTable(T emptyVal, size_t initSize) :
		valueFunc([](index) -> T {throw std::runtime_error("Set a function to memoize first!");}),
		valueFuncSet(false), emptyVal(emptyVal) {
	data.resize(initSize, emptyVal);
}

template<typename T>
bool MemoizationTable<T>::trySetValueFunc(std::function<T(index)> function) {
	if (valueFuncSet)
		return false;
	MemoizationTable::valueFunc = function;
	valueFuncSet = true;
	return true;
}

template<typename T>
T MemoizationTable<T>::getValue(index i) {
	if (data.size() <= i)
		data.resize(i + 1, emptyVal);
	if (data[i] == emptyVal)
		data[i] = valueFunc(i);
	return data[i];
}

} // namespace NetworKit

#endif //NETWORKIT_MEMOIZATIONTABLE_H
