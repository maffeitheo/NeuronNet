#include "random.h"

RandomNumbers::RandomNumbers(unsigned long int s) {
	if (seed == 0) {
		std::random_device rd;
		seed = rd();
	}
	
	rng = std::mt19937(seed);
}

void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper) {
	for (int i(0); i < vec.size(); i++) {
		vec[i] = uniform_double(lower, upper);
	}
}

double RandomNumbers::uniform_double(double lower, double upper) {
	std::uniform_real_distribution<> unif(lower, upper);
	return unif(rng);
}

void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd) {
	for (int i(0); i < vec.size(); i++) {
		vec[i] = normal(mean, sd);
	}
}

double RandomNumbers::normal(double mean, double sd) {
	std::normal_distribution<> norm(mean, sd);
	return norm(rng);
}

void RandomNumbers::poisson(std::vector<int>& vec, double mean) {
	for (int i(0); i < vec.size(); i++) {
		vec[i] = poisson(mean);
	}
}

int RandomNumbers::poisson(double mean) {
	std::poisson_distribution<> poisson(mean);
	return poisson(rng);
}

