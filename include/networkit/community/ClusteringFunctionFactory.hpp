/*
 * ClusteringFunctionFactory.h
 *
 * Created on: 2019-11-12
 * Author: Armin Wiebigke
  */

#ifndef NETWORKIT_CLUSTERINGFUNCTIONFACTORY_HPP
#define NETWORKIT_CLUSTERINGFUNCTIONFACTORY_HPP

namespace NetworKit {

using ClusteringFunction = std::function<Partition(const Graph &)>;

class ClusteringFunctionFactory {
public:
	virtual ClusteringFunction getFunction() const = 0;
};

} // namespace NetworKit

#endif //NETWORKIT_CLUSTERINGFUNCTIONFACTORY_HPP
