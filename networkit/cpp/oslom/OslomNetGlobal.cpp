#include "OslomNetGlobal.h"

OslomNetGlobal::OslomNetGlobal(std::map<int, std::map<int, std::pair<int, double> > > &A)
		: oslom_net_evaluate(A) {
}

OslomNetGlobal::OslomNetGlobal(IntMatrix &b,
                               std::deque<std::deque<std::pair<int, double> > > &c,
                               std::deque<int> &d) : oslom_net_evaluate(b, c, d) {
}

OslomNetGlobal::OslomNetGlobal(std::string a) : oslom_net_evaluate(a) {
}

/**
 * Try to merge insignificant modules to create significant modules.
 * @param discarded (in): The insignificant modules
 * @param good_modules (inout): Append the new merged modules
 * @param bscores_good (inout): Append the B-Scores of the new modules
 * @param new_discarded (out): The remaining insignificant modules
 * @return
 */
int OslomNetGlobal::try_to_merge_discarded(IntMatrix &discarded,
                                           IntMatrix &good_modules,
                                           std::deque<double> &bscores_good,
                                           IntMatrix &new_discarded) {
	// TODO: How does this work?
	// TODO: Keep merged groups even if they get extended/shrunk by more than factor 2 (?)
	new_discarded.clear();
	if (discarded.empty())
		return -1;
	if (paras->print_flag_subgraph)
		std::cout << "checking unions of not significant modules, modules to check: "
		          << discarded.size() << std::endl;

	ModuleCollection discarded_modules(dim);
	for (auto &i : discarded)
		discarded_modules.insert(i, 1);

	// Cut between discarded modules (?)
	std::map < int, std::map < int, std::pair < int, double > > >
	                                                 neigh_weight_s; // this maps the module id into the neighbor module ids and weights
	set_upper_network(neigh_weight_s, discarded_modules);

	if (neigh_weight_s.empty())
		return -1;

	OslomNetGlobal community_net(neigh_weight_s); // Graph with discarded modules as nodes
	IntMatrix M_raw;  /* M_raw contains the module_id of discarded_modules */
	community_net.collect_raw_groups_once(M_raw);

	for (auto &merged_modules : M_raw) {
		if (merged_modules.size() > 1) {
			std::set<int> l1;
			for (int module_id : merged_modules)
				deque_numeric::deque_to_set_app(discarded_modules.modules[module_id], l1);

			std::deque<int> merged_group;
			deque_numeric::set_to_deque(l1, merged_group);
//            print_id(merged_group, std::cout);
			std::deque<int> cleaned_group;
			double b_score = CUP_check(merged_group, cleaned_group);

			bool added = try_add_good_group(cleaned_group, b_score, merged_group, good_modules,
			                                bscores_good, new_discarded);
//			if (added) {
//				std::cout << "merged " << merged_modules.size() << " discarded modules "
//				          << std::endl;
//				std::cout << merged_group.size() << " -> " << cleaned_group.size() << std::endl;
//			}
		}
	}
	return 0;
}

void OslomNetGlobal::get_single_trial_partition(IntMatrix &good_modules_to_prune,
                                                std::deque<double> &bscores_good) {
	/* this function collects significant modules in two steps:
	   1. using collect_raw_groups_once and cleaning
	   2. putting together discarded modules
	   the results are appended in the input data	*/
	IntMatrix discarded;
	IntMatrix M;

	//*****************************************************************************
	collect_raw_groups_once(M);

	//cout<<"single_gather"<<endl;
	//print_id(M, std::cout);

	//*****************************************************************************

	count total_nodes_tested = 0;

	for (count i = 0; i < M.size(); i++) {

		if (paras->print_flag_subgraph && i % 100 == 0) {
			std::cout << "checked " << i << " modules " << good_modules_to_prune.size()
			          << " were found significant.  Modules to check: " << M.size() - i
			          << ". Percentage nodes done: " << double(total_nodes_tested) / dim
			          << std::endl;
		}

		/*if(paras->print_flag_subgraph && M[i].size()>1000)
			std::cout<<"M[i].size(): "<<M[i].size()<<endl;*/

		total_nodes_tested += M[i].size();

		std::deque<int> l;
		double bscore;

		if (M[i].size() < 1000)
			bscore = group_inflation(M[i], l);
		else        /* if M[i] is big enough the check is faster to save time */
			bscore = CUP_both(M[i], l);

		if (!l.empty()) {
			good_modules_to_prune.push_back(l);
			bscores_good.push_back(bscore);
		} else
			discarded.push_back(M[i]);
	}

	if (paras->print_flag_subgraph)
		std::cout << "significance check done " << std::endl << std::endl << std::endl;

	//*****************************************************************************

	IntMatrix new_discarded;
	try_to_merge_discarded(discarded, good_modules_to_prune, bscores_good, new_discarded);
	discarded = new_discarded;
	try_to_merge_discarded(discarded, good_modules_to_prune, bscores_good, new_discarded);

	if (paras->print_flag_subgraph)
		std::cout << "checking unions of not significant modules done " << std::endl << std::endl
		          << std::endl;

	/* actually here it is also possible to check the new_discarded modules more than twice but I believe this should be enough */
	/* in principle one could do while(try_to_merge_discarded(...)!=-1) */
}

/**
 *
 * @param good_modules_to_prune (out): The modules
 * @param bscores_good
 * @param runs
 */
void OslomNetGlobal::single_gather(IntMatrix &good_modules_to_prune,
                                   std::deque<double> &bscores_good, int runs = 1) {
	good_modules_to_prune.clear();
	bscores_good.clear();

	for (int i = 0; i < runs; i++)
		get_single_trial_partition(good_modules_to_prune, bscores_good);
}

/* Minimality means that the cluster has no internal structures
 * @input A : List of clusters, each cluster a list of nodes (list == deque)
 *            The input cover. (Vielleicht auch was anderes)
 */
void OslomNetGlobal::check_minimality_all(IntMatrix &A, std::deque<double> &bss,
                                          ModuleCollection &minimal_modules) {
	if (!paras->check_minimality)
		throw std::runtime_error("Check minimality is disabled!");
	paras->print_flag_subgraph = false;
	{
		/* simplifying A*/
		ModuleCollection suggestion_mall(dim);
		for (auto &i : A)
			suggestion_mall.insert(i, 1);

		suggestion_mall.erase_included();
		suggestion_mall.set_partition(A);
	}

	int counter = 0;
	while (!A.empty()) {
		IntMatrix suggestion_matrix;
		std::deque<double> suggestion_bs;

		check_minimality_matrix(A, bss, minimal_modules, suggestion_matrix, suggestion_bs,
		                        counter);

		ModuleCollection suggestion_mall(dim);
		for (auto &i : suggestion_matrix)
			suggestion_mall.insert(i, 1);

		suggestion_mall.erase_included();
		suggestion_mall.set_partition(A);

		bss = suggestion_bs;
		++counter;
	}
}

void OslomNetGlobal::check_minimality_matrix(IntMatrix &A, std::deque<double> &bss,
                                             ModuleCollection &minimal_modules,
                                             IntMatrix &suggestion_matrix,
                                             std::deque<double> &suggestion_bs,
                                             int counter) {
	if (A.size() > 4)
		std::cout << "minimality check: " << A.size() << " modules to check, run: " << counter
		          << std::endl;
	if (counter < paras->minimality_stopper) {
		for (count i = 0; i < A.size(); i++) {
			check_minimality(A[i], bss[i], minimal_modules, suggestion_matrix,
			                 suggestion_bs);
		}
	} else {
		for (count i = 0; i < A.size(); i++)
			minimal_modules.insert(A[i], bss[i]);
	}
}

/**
 * This function checks the minimality of group.
 * Minimality means that group doesn't have internal structures up to
 *   a factor coverage_percentage_fusion_left.
 * @param group (in): Module to check
 * @param bs_group (out): B-Score of the module
 * @param minimal_modules (in-out): Contains significant minimal modules. Append the input
 *                                  module if it is minimal.
 * @param suggestion_matrix (out): ?
 * @param suggestion_bs (out): The B-Scores of the suggestion_matrix modules
 * @return true iff group is inserted in minimal_modules
 */
bool OslomNetGlobal::check_minimality(std::deque<int> &group, double &bs_group,
                                      ModuleCollection &minimal_modules,
                                      IntMatrix &suggestion_matrix,
                                      std::deque<double> &suggestion_bs) {
	IntMatrix subM;
	std::deque<double> bss;
	{    //******************  module_subgraph stuff   ******************

		std::deque<std::deque<int> > link_per_node;
		std::deque<std::deque<std::pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		OslomNetGlobal module_subgraph(link_per_node, weights_per_node, group);
		std::deque<double> bscores_good_temp;
		module_subgraph.single_gather(subM, bscores_good_temp);

		for (auto &i : subM) {
			module_subgraph.deque_id(i);
			std::deque<int> grbe;
			bss.push_back(CUP_check(i, grbe));
			i = grbe;
		} /* so now you know these modules are cleaned (but you are not sure they are minimal) */
	}   //******************  module_subgraph stuff   ******************

	for (auto &i : subM)
		if (i.size() == group.size()) {
			minimal_modules.insert(group, bs_group);
			return true;
		}

	std::set<int> a;
	for (auto &i : subM)
		for (int j : i)
			a.insert(j);

	if (a.size() > paras->coverage_percentage_fusion_left * group.size()) {
		/* this means the group cannot be accepted */
		for (count i = 0; i < subM.size(); i++)
			if (!subM[i].empty()) {
				suggestion_matrix.push_back(subM[i]);
				suggestion_bs.push_back(bss[i]);
			}
		return false;
	} else {
		minimal_modules.insert(group, bs_group);
		return true;
	}
}

void
OslomNetGlobal::print_modules(bool not_homeless, const std::string &tp, ModuleCollection &Mcoll) {
	char b[1000];
	cast_string_to_char(tp, b);
	std::ofstream out1(b);
	print_modules(not_homeless, out1, Mcoll);
}

void OslomNetGlobal::print_modules(bool not_homeless, std::ostream &out1,
                                   ModuleCollection &Mcoll) {
	int nmod = 0;
	for (auto itm = Mcoll.module_bs.begin();
	     itm != Mcoll.module_bs.end(); itm++)
		if (Mcoll.modules[itm->first].size() > 1)
			nmod++;

	std::cout << "******** ModuleCollection ******** " << nmod << " modules. writing... "
	          << std::endl;

	std::deque<int> netlabs;
	for (int i = 0; i < dim; i++)
		netlabs.push_back(id_of(i));

	Mcoll.print(out1, netlabs, not_homeless);
	std::cout << "DONE   ****************************" << std::endl;
}

void OslomNetGlobal::print_statistics(std::ostream &outt, ModuleCollection &Mcoll) {
	int nmod = 0;
	count cov = 0;

	for (auto itm = Mcoll.module_bs.begin();
	     itm != Mcoll.module_bs.end(); itm++)
		if (Mcoll.modules[itm->first].size() > 1) {
			nmod++;
			cov += Mcoll.modules[itm->first].size();
		}

	std::deque<int> homel;
	Mcoll.homeless(homel);

	outt << "number of modules: " << nmod << std::endl;
	outt << "number of covered nodes: " << dim - homel.size()
	     << " fraction of homeless nodes: " << double(homel.size()) / dim << std::endl;
	outt << "average number of memberships of covered nodes: "
	     << double(cov) / (dim - homel.size()) << std::endl;
	outt << "average community size: " << double(cov) / nmod << std::endl;

	print_degree_of_homeless(homel, outt);
}

void
from_IntMatrix_and_deque_to_deque(IntMatrix &its_submodules, const std::deque<int> &A,
                                  std::deque<int> &group) {
	// it merges A and its_submodules in group
	std::set<int> all_the_groups;
	for (auto &its_submodule : its_submodules) {
		for (int j : its_submodule)
			all_the_groups.insert(j);
	}

	for (int i : A)
		all_the_groups.insert(i);

	deque_numeric::set_to_deque(all_the_groups, group);
}

bool OslomNetGlobal::fusion_module_its_subs(const std::deque<int> &A,
                                            std::deque<std::deque<int>> &its_submodules) {
	// A is supposed to be a good cluster
	// return true if A won against its submodules
	// ******************************************
	if (its_submodules.size() < 2)
		return true;

	std::deque<int> group;
	from_IntMatrix_and_deque_to_deque(its_submodules, A, group);
	{    //******************  sub_graph_module stuff   ******************
		std::deque<std::deque<int>> link_per_node;
		std::deque<std::deque<std::pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		OslomNetGlobal sub_graph_module(link_per_node, weights_per_node, group);

		sub_graph_module.translate(its_submodules);

		//------------------------------------ cleaning up submodules --------------------------
		ModuleCollection sub_mall(sub_graph_module.dim);

		for (auto &its_submodule : its_submodules)
			sub_mall.insert(its_submodule, 1e-3);

		sub_mall.set_partition(its_submodules);

		/*
		std::cout<<"group*************************************************"<<endl;
		print_id(group, std::cout);
		std::cout<<"A"<<endl;
		print_id(A, std::cout);
		std::cout<<"fusion_module_its_subs"<<endl;
		print_id(its_submodules, std::cout);
		//*/

		//------------------------------------ cleaning up submodules --------------------------
		std::set<int> a;

		for (const auto &its_submodule : its_submodules) {

			std::deque<int> grbe;
			sub_graph_module.CUP_check(its_submodule, grbe);
			deque_numeric::deque_to_set_app(grbe, a);
			//cout<<i<<" cleaned_up: "<<grbe.size()<<" "<<a.size()<<endl;

			if (a.size() > paras->coverage_percentage_fusion_or_submodules * A.size())
				return false;
		}
		//sub_graph_module.draw("sub");
		//cherr();
		return true;
	}   //******************  sub_graph_module stuff   ******************
}

bool
OslomNetGlobal::fusion_with_empty_A(IntMatrix &its_submodules, std::deque<int> &A, double &bs) {
	/*
		its_submodules are the modules to check. the question is if to take its_submodules or the union of them
		the function returns true if it's the union, grc1 is the union cleaned and bs the score
	 */
	std::deque<int> group;
	from_IntMatrix_and_deque_to_deque(its_submodules, A, group);

	//cout<<"trying a module of "<<group.size()<<" nodes"<<endl;

	bs = CUP_check(group, A);

	if (A.size() <= paras->coverage_percentage_fusion_left * group.size()) {
		A.clear();
		bs = 1;
		return false;
	}
	bool fus = fusion_module_its_subs(A, its_submodules);
	return fus;
}

void OslomNetGlobal::check_existing_unions(ModuleCollection &mall) {
	/* this function is to check unions of existing modules*/
	/* sorting from the biggest to the smallest module */
	/*cout<<"before check_existing_unions"<<endl;
	print_modules(false, std::cout, mall);*/

	std::deque<int> sm;
	mall.sort_modules(sm);

	/*cout<<"sm"<<endl;
	prints(sm);*/


	std::deque<bool> still_good;
	for (count i = 0; i < sm.size(); i++)
		still_good.push_back(true);

	std::set<int> modules_to_erase;
	for (int i : sm) {
		/* for each module I check if it's better to take it or its submodules */
		std::deque<int> smaller;
		mall.almost_equal(i, smaller);
		IntMatrix its_submodules;
		for (int j : smaller)
			if (still_good[j])
				its_submodules.push_back(mall.modules[j]);
		/*cout<<"************************** module to check "<<sm[i]<<" size: "<<mall.modules[sm[i]].size()<<endl;
		print_id(mall.modules[sm[i]], std::cout);
		std::cout<<"its_submodules"<<endl;
		print_id(its_submodules, std::cout);*/

		if (fusion_module_its_subs(mall.modules[i], its_submodules)) {
			deque_numeric::deque_to_set_app(smaller, modules_to_erase);
			for (int j : smaller)
				still_good[j] = false;
		} else {
			modules_to_erase.insert(i);
			still_good[i] = false;
		}
	}

	for (int its : modules_to_erase)
		mall.erase(its);

	mall.compact();
	/*cout<<"after check_existing_unions --------------------------------------------------------"<<endl;
	print_modules(false, std::cout, mall);*/
}

bool OslomNetGlobal::check_fusion_with_gather(ModuleCollection &mall) {
	/*	this function is used to check if we would like unions of modules
		returns true if it merges something	 */
	std::cout << "check unions of modules using community network" << std::endl << std::endl;
	paras->print_flag_subgraph = true;

	mall.fill_gaps();
	std::map<int, std::map<int, std::pair<int, double> > >
			neigh_weight_s;        // this maps the module id into the neighbor module ids and weights
	set_upper_network(neigh_weight_s, mall);

	if (neigh_weight_s.empty())
		return false;

	for (auto &neigh_weight_ : neigh_weight_s) {
		for (auto &itm2 : neigh_weight_.second) {
			itm2.second.first = 1;
		}
	}

	bool real_paras_weighted = paras->weighted;
	paras->weighted = true;
	OslomNetGlobal community_net(neigh_weight_s);
	IntMatrix M_raw;        /* M_raw contains the module_ids */
	community_net.collect_raw_groups_once(M_raw);

	paras->weighted = real_paras_weighted;
	bool something = false;
	IntMatrix module_to_insert;
	std::deque<double> bs_to_insert;
	int fused_modules = 0;
	std::cout << "possible fusions to check: " << M_raw.size() << std::endl;

	for (count i = 0; i < M_raw.size(); i++)
		if (M_raw[i].size() > 1) {
			IntMatrix ten;
			for (count j = 0; j < M_raw[i].size(); j++)
				ten.push_back(mall.modules[M_raw[i][j]]);
			//cout<<"trying fusion # "<<i<<" "<<ten.size()<<" modules to merge"<<endl;

			std::deque<int> grc1;
			double bs;
			if (fusion_with_empty_A(ten, grc1, bs)) {
				something = true;
				module_to_insert.push_back(grc1);
				++fused_modules;
				bs_to_insert.push_back(bs);
			}
//			if (i % 100 == 0)
//				std::cout << "checked " << i << " unions. Fused: " << fused_modules << std::endl;
		}

	for (count i = 0; i < module_to_insert.size(); i++)
		mall.insert(module_to_insert[i], bs_to_insert[i]);

	mall.compute_inclusions();
	return something;
}

int OslomNetGlobal::check_unions_and_overlap(ModuleCollection &mall, bool only_similar) {
	mall.put_gaps();
	if (mall.effective_groups() == 0)
		return 0;

	std::cout << "checking similar modules" << std::endl << std::endl;
	check_existing_unions(mall);
	if (!only_similar) {
		if (check_fusion_with_gather(mall))
			check_fusion_with_gather(mall);
	}
	std::cout << "checking highly intersecting modules" << std::endl << std::endl;
	check_intersection(mall);
	mall.compute_inclusions();
	return 0;
}

int OslomNetGlobal::try_to_assign_homeless(ModuleCollection &Mcoll, bool anyway) {
	Mcoll.put_gaps();
	//if(paras->print_cbs)
	//cout<<"checking homeless nodes "<<endl;

	std::deque<int> homel;
	Mcoll.homeless(homel);

	int before_procedure = homel.size();
	if (homel.empty())
		return before_procedure;

	/*cout<<"homel"<<endl;
	print_id(homel, cout);*/
	std::set<int> called;                   // modules connected to homeless nodes
	std::map<int, std::set<int> > homel_module;  // maps the homeless node with the modules it's connected to

	for (int i : homel) {
		std::set<int> thish;
		for (int j = 0; j < vertices[i]->links->size(); j++) {
			int &neigh = vertices[i]->links->l[j];
			for (int itk : Mcoll.memberships[neigh]) {
				called.insert(itk);
				thish.insert(itk);
			}
		}
		if (!thish.empty())
			homel_module[i] = thish;
	}

	std::map<int, int> module_kin;
	std::map<int, int> module_ktot;
	for (int its : called) {
		module_kin[its] = cast_int(kin_m(Mcoll.modules[its]));
		module_ktot[its] = cast_int(ktot_m(Mcoll.modules[its]));
	}

	std::map<int, std::deque<int> > to_check;            // module - homeless nodes added to that
	for (auto &itm : homel_module) {
		double cmin = 1.1;
		int belongs_to = -1;
		//cout<<"homeless node: "<<id_of(itm->first)<<endl;
		for (auto its = itm.second.begin();
		     its != itm.second.end(); its++) {
			int kin_node = cast_int(vertices[itm.first]->kplus_m(Mcoll.modules[*its]));
			/*cout<<"module: "<<*its<<" kin: "<<module_kin[*its]<<"  ktot: "<<module_ktot[*its]<<" kin h "<<kin_node<<endl;
			print_ri(Mcoll.modules[*its]);*/
			int kout_g = module_ktot[*its] - module_kin[*its];
			int tm = total_stubs - module_ktot[*its];
			//double rh= compute_r_hyper(kin_node, kout_g, tm, vertices[itm->first]->stub_number);
			double kinw = vertices[itm.first]->kplus_w(Mcoll.modules[*its]);
			//double weight_part= log_together(kinw, kin_node);
			double rh = Stochastics::compute_global_fitness_randomized_short(
					kin_node, kout_g, tm, vertices[itm.first]->stub_number, kinw);

			//double cs=  1 - pow(1 - rh, dim - Mcoll.modules[*its].size());
			//cout<<"rh: "<<rh<<" ..."<<endl;
			if (rh < cmin) {
				cmin = rh;
				belongs_to = *its;
			}
		}

		if (belongs_to != -1) {
			if (to_check.find(belongs_to) == to_check.end()) {
				std::deque<int> void_d;
				to_check[belongs_to] = void_d;
			}
			to_check[belongs_to].push_back(itm.first);
		}
		//if(paras->print_cbs)
		//cout<<"homeless node: "<<id_of(itm->first)<<" belongs_to "<<belongs_to<<" cmin... "<<cmin<<endl;
		//cherr();
	}
	//if(paras->print_cbs)
	//cout<<"homeless node: "<<homel.size()<<" try_to_assign: "<<homel_module.size()<<" modules to check: "<<to_check.size()<<endl;

	// **** try the groups with the homeless //******************
	bool something = false;
	for (auto &itm : to_check) {
		std::deque<int> union_deque = Mcoll.modules[itm.first];

		for (int i : itm.second)
			union_deque.push_back(i);
		if (anyway) {
			something = true;
			Mcoll.insert(union_deque, ran4() + paras->threshold);
		} else {
			std::deque<int> grbe;
			double bs = CUP_check(union_deque, grbe);
			//cout<<"union_deque after "<<itm->first<<" size: "<<grbe.size()<<endl;
			if (grbe.size() > 1) {
				something = true;
				Mcoll.insert(grbe, bs);
			}
		}
	}

	if (something) {
		Mcoll.compute_inclusions();
	}

	return before_procedure;
}

int OslomNetGlobal::check_intersection(ModuleCollection &Mcoll) {

	paras->print_flag_subgraph = false;

	std::deque<int> to_check;
	for (auto &module_b : Mcoll.module_bs)
		to_check.push_back(module_b.first);

	return check_intersection(to_check, Mcoll);
}

int OslomNetGlobal::check_intersection(std::deque<int> &to_check, ModuleCollection &Mcoll) {

	std::set<std::pair<int, int>> pairs_to_check;

	for (int &itM : to_check)
		if (Mcoll.module_bs.find(itM) != Mcoll.module_bs.end()) {

			std::deque<int> &c = Mcoll.modules[itM];
			std::map<int, int> com_ol;                        // it maps the index of the modules into the overlap (overlap=number of overlapping nodes)

			for (int i : c)
				for (int itj : Mcoll.memberships[i])
					int_histogram(itj, com_ol);

			for (auto &cit : com_ol)
				if (cit.first != itM) {
					if (double(cit.second) /
					    std::min(Mcoll.modules[cit.first].size(), c.size()) >
					    paras->check_inter_p) {        // they have a few nodes in common
						pairs_to_check.insert(
								std::make_pair(std::min(itM, cit.first), std::max(itM, cit.first)));
					}
				}
		}
	return fusion_intersection(pairs_to_check, Mcoll);
}

int OslomNetGlobal::fusion_intersection(std::set<std::pair<int, int>> &pairs_to_check,
                                        ModuleCollection &Mcoll) {
	std::cout << "pairs to check: " << pairs_to_check.size() << std::endl;
	std::deque<int> new_insertions;
	for (const auto &ith : pairs_to_check)
		if (ith.first < ith.second)
			if (Mcoll.module_bs.find(ith.first) != Mcoll.module_bs.end())
				if (Mcoll.module_bs.find(ith.second) != Mcoll.module_bs.end()) {
					//		first, you need to check if both the modules in the pair are still in mcoll
					std::deque<int> &a1 = Mcoll.modules[ith.first];
					std::deque<int> &a2 = Mcoll.modules[ith.second];
					int min_s = std::min(a1.size(), a2.size());
					std::deque<int> group_intsec;
					std::set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(),
					                      back_inserter(group_intsec));
					//		if they are, you need to check if they are not almost equal.
					if (double(group_intsec.size()) / min_s >=
					    paras->coverage_inclusion_module_collection) {
						int em = ith.first;
						if (a1.size() < a2.size())
							em = ith.second;
						else if (a1.size() == a2.size() &&
						         Mcoll.module_bs[ith.first] > Mcoll.module_bs[ith.second])
							em = ith.second;
						Mcoll.erase(em);
					} else
						decision_fusion_intersection(ith.first, ith.second, new_insertions,
						                             Mcoll,
						                             double(group_intsec.size()) / min_s);
				}
	if (!new_insertions.empty())
		return check_intersection(new_insertions, Mcoll);
	return 0;
}

bool OslomNetGlobal::decision_fusion_intersection(int ai1, int ai2,
                                                  std::deque<int> &new_insertions,
                                                  ModuleCollection &Mcoll,
                                                  double prev_over_percentage) {
	std::deque<int> &a1 = Mcoll.modules[ai1];
	std::deque<int> &a2 = Mcoll.modules[ai2];

	std::deque<int> group;
	std::set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group));

	if (int(group.size()) != dim) {
		//******************  sub_graph_module stuff   ******************

		std::deque<std::deque<int>> link_per_node;
		std::deque<std::deque<std::pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		OslomNetGlobal sub_graph_module(link_per_node, weights_per_node, group);

		std::deque<std::deque<int>> A;
		A.push_back(a1);
		A.push_back(a2);

		sub_graph_module.translate(A);

		std::deque<int> grc1;
		double bs = sub_graph_module.CUP_check(A[0], grc1);
		std::deque<int> grc2;
		bs = sub_graph_module.CUP_check(A[1], grc2);

		std::deque<int> unions_grcs;
		std::set_union(grc1.begin(), grc1.end(), grc2.begin(), grc2.end(),
		               back_inserter(unions_grcs));

		if (unions_grcs.size() <= paras->coverage_percentage_fusion_or_submodules *
		                          group.size()) {        // in such a case you can take the fusion (if it's good)

			/* actually the right check should be  unions_grcs.size() > paras->coverage_percentage_fusion_or_submodules*group_2.size()
			   but this would require more time - it should not make a big difference anyway */
			std::deque<int> group_2;
			bs = CUP_check(group, group_2);
			if (group_2.size() > paras->coverage_percentage_fusion_left * group.size()) {
				Mcoll.erase(ai1);
				Mcoll.erase(ai2);
				IntMatrix _A_;
				std::deque<double> _bss_;
				_A_.push_back(group_2);
				_bss_.push_back(bs);
				check_minimality_all(_A_, _bss_, Mcoll);
				return true;
			} else
				return false;
		}

		sub_graph_module.deque_id(grc1);
		sub_graph_module.deque_id(grc2);
		std::deque<int> cg1;
		double bs__1 = CUP_check(grc1, cg1);
		std::deque<int> cg2;
		double bs__2 = CUP_check(grc2, cg2);
		std::deque<int> inters;
		std::set_intersection(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(),
		                      back_inserter(inters));
		std::deque<int> unions;
		std::set_union(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(), back_inserter(unions));

		if (double(inters.size()) / std::min(cg1.size(), cg2.size()) <
		    prev_over_percentage - 1e-4) {
			if (!cg1.empty() && !cg2.empty() &&
			    (unions.size() > paras->coverage_percentage_fusion_left * group.size())) {
				Mcoll.erase(ai1);
				Mcoll.erase(ai2);
				int newi;
				Mcoll.insert(cg1, bs__1, newi);
				new_insertions.push_back(newi);
				Mcoll.insert(cg2, bs__2, newi);
				new_insertions.push_back(newi);
				//cout<<"pruned module"<<endl;
				return true;
			}
		}
	}
	return false;
}

void print_seperator_line() {
	std::cout << "***************************************************************************"
	          << std::endl;
}

/**
 * Try to add a cleaned group to the good modules. If the group size increased or shrunk too much,
 * discard the group.
 * @param group
 * @param b_score
 * @param original_group
 * @param good_modules
 * @param bscores_good
 * @param bad_groups
 * @return
 */
bool OslomNetGlobal::try_add_good_group(
		std::deque<int> &group, double &b_score, const std::deque<int> &original_group,
		IntMatrix &good_modules, std::deque<double> &bscores_good, IntMatrix &bad_groups) {
	// Just remove bad nodes if the group is extended too much
	if (!group.empty() && group.size() > paras->max_group_extend * original_group.size()) {
		group.clear();
		if (!paras->discard_max_extend_groups) {
			throw std::runtime_error("Discard");
			b_score = CUP_check(original_group, group);
//			std::cout << original_group.size() << " -> " << group.size() << std::endl;
		}
	}
	bool added = true;
	// A group is only good if it wasn't shrunk too much
	if (!group.empty() && group.size() >= original_group.size() / paras->max_group_extend) {
		good_modules.push_back(group);
		bscores_good.push_back(b_score);
	} else {
//		std::cout << "=> bad group" << std::endl;
		if (paras->keep_bad_groups) {
			good_modules.push_back(original_group);
			bscores_good.push_back(0);
		} else {
			bad_groups.push_back(original_group);
			added = false;
		}
	}
	return added;
}

ModuleCollection
OslomNetGlobal::clean_up(const std::vector<std::deque<int>> &modules, int upper_node_id) {
	ModuleCollection minimal_modules(upper_node_id);
	IntMatrix good_modules;
	std::deque<double> bscores_good;
	std::ofstream bad_groups_file(paras->bad_groups_filename);
	IntMatrix bad_groups;

	std::cout << modules.size() << " groups found" << std::endl;
	for (const auto & module : modules) {
		std::deque<int> group;
//        std::cout << "processing group number " << gr_id << " size: " << modules[gr_id].size()
//                  << std::endl;
		double bcu = 0.0;
		if (paras->cleanup_strategy == "both")
			bcu = CUP_both(module, group);
		else if (paras->cleanup_strategy == "check")
			bcu = CUP_check(module, group);
		else if (paras->cleanup_strategy == "search")
			bcu = CUP_search(module, group);
		else
			throw std::runtime_error("No cleanup strategy!");
//        std::cout << modules[gr_id].size() << " -> " << group.size() << std::endl;

		try_add_good_group(group, bcu, module, good_modules, bscores_good, bad_groups);
	}
	print_seperator_line();


	for (const auto &bad_group : bad_groups) {
		for (int u : bad_group)
			bad_groups_file << u << " ";
		bad_groups_file << std::endl;
	}

	// Merge discarded groups
	if (paras->merge_discarded) {
		std::cout << "TRY TO MERGE DISCARDED GROUPS" << std::endl;
		IntMatrix remaining_bad_groups;
		try_to_merge_discarded(bad_groups, good_modules, bscores_good,
		                       remaining_bad_groups);
		bad_groups = remaining_bad_groups;
		try_to_merge_discarded(bad_groups, good_modules, bscores_good,
		                       remaining_bad_groups);
		std::cout << "MERGE DISCARED GROUPS DONE" << std::endl;
		print_seperator_line();
	}

	if (paras->check_minimality) {
		// Ensure that modules are minimal
		check_minimality_all(good_modules, bscores_good, minimal_modules);
		std::cout << "***************************************************************************"
		          << std::endl;
		std::cout << "MINIMALITY CHECK DONE" << std::endl;
	} else {
		for (size_t i = 0; i < good_modules.size(); ++i) {
			minimal_modules.insert(good_modules[i], bscores_good[i]);
		}
	}

	// Try to merge modules
	if (paras->check_unions) {
		check_unions_and_overlap(minimal_modules);
		std::cout << "***************************************************************************"
		          << std::endl;
		std::cout << "CHECK UNIONS AND SIMILAR MODULES DONE" << std::endl;
	}
	return minimal_modules;
}


