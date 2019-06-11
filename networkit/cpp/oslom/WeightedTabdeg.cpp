#include "WeightedTabdeg.h"
#include "Stochastics.h"

void WeightedTabdeg::clear() {
	labels_to_facts.clear();
	fitness_to_label.clear();
}

/**
 * Insert an element. If it was already inserted, edit the element instead.
 * @param a Node label
 * @param kp Internal Degree
 * @param kt Degree
 * @param mtlw minus_log_total_wr
 * @param fit Fitness Score
 */
void WeightedTabdeg::insert_node(int a, int kp, int kt, double mtlw, double fit) {
	// this function inserts element a (or edit it if it was already inserted)
	erase_node(a);
	auto fiit = fitness_to_label.insert(std::make_pair(fit, a));
	facts F(kp, mtlw, fiit, kt);
	labels_to_facts.insert(std::make_pair(a, F));
}

/**
 * Erase an element if it exists.
 * @param a The element to earase
 * @return True iff the element exists
 */
bool WeightedTabdeg::erase_node(int a) {
	auto itm = labels_to_facts.find(a);
	if (itm != labels_to_facts.end()) {
		fitness_to_label.erase(itm->second.fitness_iterator);
		labels_to_facts.erase(itm);
		return true;
	}
	return false;
}

bool WeightedTabdeg::is_internal(int a) {
	auto itm = labels_to_facts.find(a);
	if (itm == labels_to_facts.end())
		return false;
	return true;
}

void WeightedTabdeg::set_deque(std::deque<int> &vv) {
	vv.clear();
	for (auto &lab_fact : labels_to_facts)
		vv.push_back(lab_fact.first);
}

void WeightedTabdeg::print_nodes(std::ostream &outb, std::deque<int> &lab_id) {
	std::cout << "printing nodes:.. (lab intk mtlw fitness degree) " << size() << std::endl;
	for (auto &lab_fact : labels_to_facts)
		std::cout << lab_id[lab_fact.first] << " " << lab_fact.second.internal_degree << " "
		          << lab_fact.second.internal_edgeweight << " "
		          << (lab_fact.second.fitness_iterator)->first << " " << lab_fact.second.degree
		          << std::endl;
}

int WeightedTabdeg::worst_node(int &lab, double &worst_fitness, int kout_g, int Nstar,
                               int nneighs, int tm) {
	//std::cout<<"worst_node fitness - lab - (cgroup)"<<endl;
	//prints(fitness_lab);
	lab = -1;
	worst_fitness = -1;
	auto bit = fitness_to_label.end();
	if (bit == fitness_to_label.begin())
		return -1;

	int stopper = 0;
	while (bit != fitness_to_label.begin()) {
		bit--;
		auto itm = labels_to_facts.find(bit->second);
		double F = Stochastics::compute_global_fitness_randomized(
				itm->second.internal_degree,
				kout_g + 2 * itm->second.internal_degree - itm->second.degree,
				tm + itm->second.degree,
				itm->second.degree,
				itm->second.internal_edgeweight,
				nneighs + 1, Nstar + 1);

		if (F > worst_fitness) {
			worst_fitness = F;
			lab = itm->first;
		}

		stopper++;
		if (stopper == num_up_to)
			break;
	}
	return 0;
}

int WeightedTabdeg::best_node(int &lab, double &best_fitness, int kout_g, int Nstar,
                              int nneighs, int tm) {
	// I can try to compute the fitness here
	/*std::cout<<"NE BEST NODE "<<endl;
	std::cout<<"fitness_lab  "<<endl;
	prints(fitness_lab);*/
	lab = -1;
	best_fitness = 1;

	auto bit = fitness_to_label.begin();
	if (bit == fitness_to_label.end()) {
		return -1;
	}

	int stopper = 0;
	while (bit != fitness_to_label.end()) {
		auto itm = labels_to_facts.find(bit->second);
		double F = Stochastics::compute_global_fitness_randomized(
				itm->second.internal_degree, kout_g,
				tm, itm->second.degree,
				itm->second.internal_edgeweight,
				nneighs, Nstar);
		//std::cout<<itm->first<<" "<<F<<" ... node-fit"<<endl;
		if (F < best_fitness) {
			best_fitness = F;
			lab = itm->first;
		}

		stopper++;
		if (stopper == num_up_to)
			break;

		bit++;
	}

	return 0;
}

void WeightedTabdeg::_set_(WeightedTabdeg &one) {
	clear();
	for (auto &lab_fact : one.labels_to_facts)
		insert_node(lab_fact.first, lab_fact.second.internal_degree, lab_fact.second.degree,
		            lab_fact.second.internal_edgeweight, (lab_fact.second.fitness_iterator)->first);
}

/**
 * Update a node after a change.
 * @param node The node to update.
 * @param delta_degree The number of edges the node gains into the group because of the change.
 * @param delta_mtlw The sum of the edgeweights the node gains into the group because of the change.
 * @param nstar
 * @param nn
 * @param kout_g
 * @param tm
 * @param kt
 * @param to_be_erased
 * @return
 */
bool WeightedTabdeg::update_group(int node, int delta_degree, double delta_mtlw, int nstar,
                                  int nn, int kout_g, int tm, int kt,
                                  std::deque<int> &to_be_erased) {
	// this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
	auto itm = labels_to_facts.find(node);
	if (itm == labels_to_facts.end())
		return false;

	itm->second.internal_edgeweight += delta_mtlw;
	itm->second.internal_degree += delta_degree;
	if (itm->second.internal_degree == 0 && size() > 1) {
		to_be_erased.push_back(node);
		return true;
	}

	//std::cout<<"UPdating... group "<<a<<endl;
	double fit = Stochastics::compute_global_fitness_ofive(
			itm->second.internal_degree,
			kout_g + 2 * itm->second.internal_degree -
			itm->second.degree,
			tm + itm->second.degree, itm->second.degree,
			itm->second.internal_edgeweight, nn + 1,
			nstar + 1);

	fitness_to_label.erase(itm->second.fitness_iterator);
	auto fiit = fitness_to_label.insert(std::make_pair(fit, node));
	itm->second.fitness_iterator = fiit;

	return true;
}

bool WeightedTabdeg::update_neighs(int a, int delta_degree, double delta_mtlw, int nstar,
                                   int kout_g, int tm, int kt) {
	// this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
	//std::cout<<"UPdating... neighs "<<a<<" "<<kt<<endl;
	auto itm = labels_to_facts.find(a);
	if (itm == labels_to_facts.end()) {
		insert_node(a, 0, kt, 0, 1);
		itm = labels_to_facts.find(a);
	}

	itm->second.internal_degree += delta_degree;
	if (itm->second.internal_degree == 0) {
		//std::cout<<"erased from neigh update "<<a<<std::endl;
		erase_node(a);
		return true;
	}

	itm->second.internal_edgeweight += delta_mtlw;
	double fit = Stochastics::compute_global_fitness_ofive(
			itm->second.internal_degree, kout_g, tm,
			itm->second.degree,
			itm->second.internal_edgeweight, size(),
			nstar);

	fitness_to_label.erase(itm->second.fitness_iterator);
	auto fiit = fitness_to_label.insert(std::make_pair(fit, a));
	itm->second.fitness_iterator = fiit;
	return true;
}

void WeightedTabdeg::set_and_update_group(int nstar, int nn, int kout_g, int tm,
                                          WeightedTabdeg &one) {
	/*this function is to set and update the fitnesses of all the nodes in cgroup*/
	clear();
	for (auto &lab_fact : one.labels_to_facts) {
		facts &node_facts = lab_fact.second;
		double fit = Stochastics::compute_global_fitness_ofive(
				node_facts.internal_degree,
				kout_g + 2 * node_facts.internal_degree - node_facts.degree,
				tm + node_facts.degree,
				node_facts.degree,
				node_facts.internal_edgeweight,
				nn + 1,
				nstar + 1);
		int node_label = lab_fact.first;
		insert_node(node_label, node_facts.internal_degree, node_facts.degree,
		            node_facts.internal_edgeweight, fit);
	}
}

void WeightedTabdeg::set_and_update_neighs(int nstar, int nn, int kout_g, int tm,
                                           WeightedTabdeg &one) {
	/*this function is to set and update the fitnesses of all the nodes in neighs*/
	clear();
	for (auto &lab_fact : one.labels_to_facts) {
		facts &node_facts = lab_fact.second;
		double fit = Stochastics::compute_global_fitness_ofive(
				node_facts.internal_degree,
				kout_g,
				tm,
				node_facts.degree,
				node_facts.internal_edgeweight,
				nn,
				nstar);
		int node_label = lab_fact.first;
		insert_node(node_label, node_facts.internal_degree, node_facts.degree,
		            node_facts.internal_edgeweight, fit);
	}
}

int WeightedTabdeg::size() { return labels_to_facts.size(); }

