//
//  particle.cpp
//

#include "../include/particle.hpp"
#include "yaml-cpp/yaml.h"
#include "../include/constants.h"
#include <vector>
#include "../exceptions.hpp"

using namespace std;
using namespace kappa;


Particle::Particle(const std::string &name, const std::string &filename) {
	readData(name, filename);
}

Particle::Particle() {};

void Particle::readData(const std::string &name, const std::string &filename) {
	YAML::Node file;

	try {
		file = YAML::LoadFile(filename);
	}
	catch (const YAML::BadFile &e) {
		std::string error_string = "Could not load database file " + filename;
		throw UnopenedFileException(error_string);
	}
    
    if (!file[name]) {
		std::string error_string = "No data found for " + name + " in the database";
		throw DataNotFoundException(error_string);
    }

    this->name = name;
    YAML::Node particle = file[name];
    if(!particle["Mass, kg"]){
        cerr << "Error: in file:" << filename << " no data for: " << name <<" : mass"<< endl;
        // FIXME
        // Need exeption!!!
        exit(-1);
    }
    // FIXME need do some if not exist in file
    if(particle["Mass, kg"]) mass =  particle["Mass, kg"].as<double>();
    if(particle["Diameter, m"]) diameter = particle["Diameter, m"].as<double>();
    if(particle["Formation energy, J"]) formation_energy = particle["Formation energy, J"].as<double>();
    if(particle["Charge"]) charge = particle["Charge"].as<int>();

    if (name != "e-") {
        if(particle["Parameter ε (Lennard-Jones), J"]) LennardJones_epsilon = particle["Parameter ε (Lennard-Jones), J"].as<double>();
        if(particle["Ionization potential, J"]) ionization_potential = particle["Ionization potential, J"].as<double>();
        if(particle["Electronic energy, J"]) electron_energy = particle["Electronic energy, J"].as<vector<double>>();
        if(particle["Statistical weight"]) statistical_weight = particle["Statistical weight"].as<vector<unsigned long>>();
        num_electron_levels = statistical_weight.size();
    }
    
}







