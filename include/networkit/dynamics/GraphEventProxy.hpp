/*
 * GraphEventProxy.h
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_GRAPH_EVENT_PROXY_HPP_
#define NETWORKIT_DYNAMICS_GRAPH_EVENT_PROXY_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/dynamics/GraphEventHandler.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 * This class enables the observer pattern for dynamic graphs: It has the same modifiers as a Graph object.
 * When these modifiers are called, they are also called on the underlying graphs. Also, all registered
 * observers (type GraphEventHandler) are notified.
 */
class GraphEventProxy {

protected:

    std::vector<GraphEventHandler*> observers;


public:

    Graph* G;

    GraphEventProxy(); // nullary constructor needed for python interface

    GraphEventProxy(Graph& G);

    void registerObserver(GraphEventHandler* observer);

    node addNode();

    void removeNode(node u);

    void restoreNode(node u);

    void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

    void removeEdge(node u, node v);

    void setWeight(node u, node v, edgeweight w);

    void incrementWeight(node u, node v, edgeweight delta);

    void TLX_DEPRECATED(timeStep());
};

} /* namespace NetworKit */
#endif // NETWORKIT_DYNAMICS_GRAPH_EVENT_PROXY_HPP_