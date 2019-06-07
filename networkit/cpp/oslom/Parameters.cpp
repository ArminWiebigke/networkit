#include <vector>

#include "Parameters.h"

void remove_element(std::vector<std::string> &arguments, int index) {
	for (size_t i = index; i < arguments.size() - 1; ++i) {
		arguments[i] = arguments[i + 1];
	}
	arguments.pop_back();
}

template<typename T>
T get_value(const std::string &value);

template<>
std::string get_value<std::string>(const std::string &value) {
	return value;
}

template<>
double get_value<double>(const std::string &value) {
	return std::stod(value);
}

template<>
int get_value<int>(const std::string &value) {
	return std::stoi(value);
}

template<typename T>
T get_arg_value(std::vector<std::string> &arguments, const std::string &name,
                T default_value) {
	size_t idx = std::find(arguments.begin(), arguments.end(), name) - arguments.begin();
	if (idx == arguments.size())
		return default_value;
	if (idx + 1 == arguments.size())
		throw std::runtime_error("Missing value for parameter \"" + name + "\"!");
	T return_val = get_value<T>(arguments[idx + 1]);
	remove_element(arguments, idx);
	remove_element(arguments, idx);
	return return_val;
}

template<>
bool get_arg_value<bool>(std::vector<std::string> &arguments, const std::string &name,
                         bool) {
	auto it = std::find(arguments.begin(), arguments.end(), name);
	if (it != arguments.end()) {
		remove_element(arguments, it - arguments.begin());
		return true;
	}
	return false;
}

bool get_bool_arg(std::vector<std::string> &arguments, const std::string &name) {
	return get_arg_value<bool>(arguments, name, false);
}

void Parameters::set(const std::vector<std::string> &args) {
	auto args_consume = args;
	check_unions = get_bool_arg(args_consume, "-check_unions");
	clean_up_runs = get_arg_value<int>(args_consume, "-cup_runs", 1);
	cleanup_strategy = get_arg_value<std::string>(args_consume, "-cu_strat", "both");
	check_minimality = get_bool_arg(args_consume, "-check_min");
	simple_cleanup = get_bool_arg(args_consume, "-simple_cleanup");
	max_group_extend = get_arg_value<double>(args_consume, "-max_extend", 1000.0);
	bad_groups_filename = get_arg_value<std::string>(args_consume, "-bad_groups_file",
	                                                 "bad_groups.txt");
	merge_discarded = get_bool_arg(args_consume, "-merge_discarded");
	threshold = get_arg_value<double>(args_consume, "-threshold", 0.1);
	keep_bad_groups = get_bool_arg(args_consume, "-keep_bad_groups");
	discard_max_extend_groups = get_bool_arg(args_consume, "-discard_max_extend_groups");
	if (!args_consume.empty()) {
		std::cout << "Remaining arguments: " << std::endl;
		for (const auto& arg : args_consume)
			std::cout << arg << " ";
		std::cout << std::endl;
		throw std::runtime_error("Unused arguments!");
	}
}

void Parameters::print() {
	std::cout << "**************************************" << std::endl;
	std::cout << "Threshold:\t\t\t" << threshold << std::endl;

	if (weighted)
		std::cout << "Weighted: yes" << std::endl;
	else
		std::cout << "Weighted: no" << std::endl;

	if (fast)
		std::cout << "-fast option selected" << std::endl;

	std::cout << "-cp:\t\t\t" << coverage_percentage_fusion_or_submodules << std::endl;

	if (!assign_homeless)
		std::cout << "-singlet option selected" << std::endl;

	std::cout << "**************************************" << std::endl << std::endl;
}

Parameters::Parameters() {
	threshold = 0.1;                                 // this is the P-value for the significance of the module

	clean_up_runs = 25;                              // the number of runs in the clean up procedure
	inflate_runs = 3;                                // the number of runs in the clean up of the inflate procedure
	inflate_stopper = 5;                             // the number of runs in the inflate procedure
	equivalence_parameter = 0.33;                    // this parameters tells when nodes are considered equivalent in the clean up procedure
	cut_off = 200;                                   // this is used in the inflate function

	maxborder_nodes = 100;                           // this is to speed up the code in looking for "reasonably good" neighbors
	maxbg_ordinary = 0.1;                            // same as above
	iterative_stopper = 10;                          // this is to prevent the iterative procedure to last too long. this can happen in case of strong backbones (just an idea, not sure)
	minimality_stopper = 10;                         // this is to prevent too many minimality checks

	coverage_inclusion_module_collection = 0.49999;  // this is used to see if two modules are higly similar in processing the clusters (big_module)
	coverage_percentage_fusion_left = 0.8;           // this is used to see when fusing clusters how much is left
	check_inter_p = 0.05;                            // this parameter is a check parameter for the fusion of quite similar clusters
	coverage_percentage_fusion_or_submodules = 0.5;  // this is the resolution parameter to decide between split clusters or unions, if you increase this value the program tends to find bigger clusters


	print_flag_subgraph = true;                      // this flag is used to print things when necessary

	fast = false;
	weighted = false;
	assign_homeless = true;

	//********************* collect_groups
	max_iteration_convergence = 10;                  // parameter for the convergence of the collect_groups function

	set({});
}

Parameters *Parameters::get_instance() {
	if (instance == nullptr)
		instance = new Parameters();
	return instance;
}

Parameters *Parameters::instance = nullptr;

void Parameters::delete_instance() {
	delete instance;
	instance = nullptr;
}

