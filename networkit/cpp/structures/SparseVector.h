/*
 * SparseVector.h
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_SPARSEVECTOR_H
#define NETWORKIT_SPARSEVECTOR_H

#include <vector>

#include "../Globals.h"

namespace NetworKit {

/**
 * A vector that imitates a map with unsigned integer keys.
 * This class has faster access than a map, but needs space linear in the maximum key value.
 */
template<typename T>
class SparseVector {
public:
	SparseVector();

	explicit SparseVector(index size);

	SparseVector(index size, T defaultValue);

	/**
	 * Resize the vector so that indexes up to size-1 can be used.
	 */
	void resize(index size);

	/**
	 * @return the size, i.e. the upper bound of the index that can be used
	 */
	index size() const;

	/**
	 * Reset all values to the default value, so it is "empty". The size is not changed.
	 */
	void reset();

	/**
	 * Returns the number of inserted elements.
	 */
	count elementCount() const;

	/**
	 * Insert an value at position a given position.
	 * @param i index where the value is inserted
	 * @param value
	 */
	void insert(index i, T value);

	/**
	 * Access operator. Before accessing an element, insert it by using the insert() method.
	 */
	T& operator[](index i);

	/**
	 * Const access operator. Before accessing an element, insert it by using the insert() method.
	 */
	const T& operator[](index i) const;

	/**
	 * Get all indexes into which elements were inserted.
	 */
	std::vector<index> insertedIndexes();

	/**
	 * Clear the value, setting the size to 0.
	 */
	void clear();

	bool indexIsUsed(index idx);

private:
	std::vector<T> data;
	std::vector<index> usedIndexes;
	T emptyValue;
};


template<typename T>
SparseVector<T>::SparseVector() : emptyValue(T{}) {
}

template<typename T>
SparseVector<T>::SparseVector(count size) : data(size), emptyValue(T{}) {
}

template<typename T>
SparseVector<T>::SparseVector(count size, T defaultValue) : data(size, defaultValue),
                                                            emptyValue(defaultValue) {
}

template<typename T>
void SparseVector<T>::reset() {
	for (index i : usedIndexes) {
		data[i] = emptyValue;
	}
	usedIndexes.clear();
}

template<typename T>
void SparseVector<T>::insert(index i, T value) {
	usedIndexes.push_back(i);
	data[i] = std::move(value);
}

template<typename T>
T& SparseVector<T>::operator[](index i) {
	return data[i];
}

template<typename T>
const T &NetworKit::SparseVector<T>::operator[](NetworKit::index i) const {
	return data[i];
}

template<typename T>
count SparseVector<T>::elementCount() const{
	return usedIndexes.size();
}

template<typename T>
void SparseVector<T>::resize(index size) {
	data.resize(size, emptyValue);
}

template<typename T>
index SparseVector<T>::size() const {
	return data.size();
}

template<typename T>
std::vector<index> NetworKit::SparseVector<T>::insertedIndexes() {
	return usedIndexes;
}

template<typename T>
void NetworKit::SparseVector<T>::clear() {
	usedIndexes.clear();
	usedIndexes.shrink_to_fit();
	data.clear();
	data.shrink_to_fit();
}

template<typename T>
bool NetworKit::SparseVector<T>::indexIsUsed(index idx) {
	return data[idx] != emptyValue;
}

} /* namespace NetworKit */

#endif //NETWORKIT_SPARSEVECTOR_H
