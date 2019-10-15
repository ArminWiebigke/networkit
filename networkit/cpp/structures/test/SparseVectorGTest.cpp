/*
 * SparseVectorGTest.cpp
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include "../SparseVector.h"


namespace NetworKit {

class SparseVectorGTest : public testing::Test {
};

TEST_F(SparseVectorGTest, testEmpty){
	SparseVector<int> vector;

	ASSERT_EQ(vector.elementCount(), 0);
}

TEST_F(SparseVectorGTest, testInsertion){
	SparseVector<int> vector(11);
	vector.insert(1, 10);
	vector.insert( 3, 11);
	vector.insert(10, 12);

	ASSERT_EQ(vector.elementCount(), 3);
	ASSERT_EQ(vector[0], 0);
	ASSERT_EQ(vector[1], 10);
	ASSERT_EQ(vector[3], 11);
	ASSERT_EQ(vector[10], 12);
}

TEST_F(SparseVectorGTest, testAssignment){
	SparseVector<int> vector(2);
	vector.insert(1, 1);
	vector[1] = 2;

	ASSERT_EQ(vector.elementCount(), 1);
	ASSERT_EQ(vector[0], 0);
	ASSERT_EQ(vector[1], 2);
}

TEST_F(SparseVectorGTest, testReset){
	SparseVector<int> vector(5);
	vector.insert(1, 1);
	vector.insert(4, 2);
	vector.reset();

	ASSERT_EQ(vector.elementCount(), 0);
	ASSERT_EQ(vector[1], 0);
	ASSERT_EQ(vector[4], 0);
}

TEST_F(SparseVectorGTest, testResize){
	SparseVector<int> vector(5);
	vector.insert(1, 1);
	vector.resize(10);

	ASSERT_EQ(vector.elementCount(), 1);
	ASSERT_EQ(vector[1], 1);
	ASSERT_EQ(vector[5], 0);
	ASSERT_EQ(vector.size(), 10);
}

TEST_F(SparseVectorGTest, testDefaultValue){
	SparseVector<int> vector(5, 1);
	vector.insert(1, 5);
	vector.resize(10);
	vector.reset();

	ASSERT_EQ(vector.elementCount(), 0);
	ASSERT_EQ(vector[1], 1);
	ASSERT_EQ(vector[5], 1);
	ASSERT_EQ(vector[9], 1);
}

}