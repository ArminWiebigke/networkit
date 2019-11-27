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

	/**
	 * Construct an empty vector. Empty values are created using the default constructor.
	 * @param size upper bound for the maximum usable index
	 */
	explicit SparseVector(index size);

	/**
	 * Construct an empty vector.
	 * @param size upper bound for the maximum usable index
	 * @param emptyValue value used for empty entries
	 */
	SparseVector(index size, T emptyValue);

	/**
	 * Resize the vector so that indexes up to size-1 can be used.
	 */
	void setUpperBound(index size);

	/**	 *
	 * @return the upper bound of the indexes that can be used
	 */
	index upperBound() const;

	/**
	 * @return the number of inserted elements.
	 */
	count size() const;

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
	const std::vector<index>& insertedIndexes() const;

	/**
	 * Returns true iff an element was previously inserted at the given index.
	 * @param idx
	 */
	bool indexIsUsed(index idx);

	/**
	 * Reset all values to the default value, so it is "empty". The upper bound is not changed.
	 */
	void reset();

	/**
	 * Clear the vector, setting the upper bound of usable indexes to 0.
	 */
	void clear();
	
	/**
	 * Reallocate the datastructure if size exceeds current upper bound
	 * This is different from setUpperBound() since we want to make sure both usedIndexes and data are allocated on the socket of the calling thread
	 * @param size
	 * @param emptyValue new emptyValue
	 */
	void resize(size_t size, T emptyValue);

private:
	std::vector<T> data;
	std::vector<index> usedIndexes;
	T emptyValue;
};


template<typename T>
SparseVector<T>::SparseVector() : emptyValue(T{}) {
}

template<typename T>
SparseVector<T>::SparseVector(count size) : SparseVector(size, T{}) {
}

template<typename T>
SparseVector<T>::SparseVector(count size, T emptyValue) : data(size, emptyValue),
                                                            emptyValue(emptyValue) {
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
	assert(i < data.size());
	return data[i];
}

template<typename T>
void SparseVector<T>::setUpperBound(index size) {
	data.resize(size, emptyValue);
}

template<typename T>
index SparseVector<T>::upperBound() const {
	return data.size();
}

template<typename T>
count SparseVector<T>::size() const{
	return usedIndexes.size();
}

template<typename T>
const std::vector<index>& NetworKit::SparseVector<T>::insertedIndexes() const {
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
void NetworKit::SparseVector<T>::resize(size_t size, T emptyValue) {
	if (size > upperBound()) {
		this->emptyValue = emptyValue;
		data = std::vector<T>(size, this->emptyValue);
		usedIndexes = std::vector<index>();
	}
}

template<typename T>
bool NetworKit::SparseVector<T>::indexIsUsed(index idx) {
	return data[idx] != emptyValue;
}

} /* namespace NetworKit */

#endif //NETWORKIT_SPARSEVECTOR_H
