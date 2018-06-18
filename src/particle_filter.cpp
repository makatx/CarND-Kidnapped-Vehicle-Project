/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <sstream>
#include <string>
#include <iterator>
#include <cmath>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;

	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);


	for (int i = 0; i < num_particles; ++i) {
		Particle p;

		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
	}
	is_initialized = true;
	std::cout << "Initialized Particle Filter..." <<endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	if(fabs(yaw_rate) > 0.001) {
		for(int i=0; i<num_particles; i++) {
			Particle p = particles[i];
			double x = p.x + velocity * (sin(p.theta+yaw_rate*delta_t) - sin(p.theta)) / yaw_rate;
			double y = p.y + velocity * (-cos(p.theta+yaw_rate*delta_t) + cos(p.theta)) / yaw_rate;
			double theta = p.theta + yaw_rate * delta_t;

			normal_distribution<double> dist_x(x, std_pos[0]);
			normal_distribution<double> dist_y(y, std_pos[1]);
			normal_distribution<double> dist_theta(theta, std_pos[2]);

			particles[i].x = dist_x(gen);
			particles[i].y = dist_y(gen);
			particles[i].theta = dist_theta(gen);
		}
	}
	else {
		for(int i=0; i<num_particles; i++) {
			Particle p = particles[i];
			double x = p.x + velocity * cos(p.theta) * delta_t;
			double y = p.y + velocity * sin(p.theta) * delta_t;

			normal_distribution<double> dist_x(x, std_pos[0]);
			normal_distribution<double> dist_y(y, std_pos[1]);

			particles[i].x = dist_x(gen);
			particles[i].y = dist_y(gen);

		}
	}
	std::cout << "Prediction complete..." <<endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//std::cout << "Number of observations: "<<observations.size()<<endl;
	//std::cout << "Number of landmarks: "<<map_landmarks.landmark_list.size()<<endl;

	for (Particle& p : particles) {
		//Transform Obeservationt from vehicle to map coordinates first:
		std::vector<LandmarkObs> obs_on_map;
		p.sense_x.clear();
		p.sense_y.clear();
		p.associations.clear();
		for(LandmarkObs l : observations) {
			LandmarkObs t;
			t.x = p.x + l.x*cos(p.theta) - l.y*sin(p.theta);
			t.y = p.y + l.x*sin(p.theta) + l.y*cos(p.theta);
			obs_on_map.push_back(t);
		}
		//std::cout << "Number of transformed observations: "<<obs_on_map.size()<<endl;

		for(LandmarkObs l : obs_on_map) {
			Map::single_landmark_s nearest_lm;
			double nearest_distance = sensor_range;
			for(Map::single_landmark_s map_lm : map_landmarks.landmark_list) {
				double distance = dist(double(map_lm.x_f), double(map_lm.y_f), l.x, l.y);
				if(distance < nearest_distance) {
					nearest_distance = distance;
					nearest_lm = map_lm;
				}
			}
			//std::cout<<"Nearest distance: "<<nearest_distance<<endl;
			p.sense_x.push_back(l.x);
			p.sense_y.push_back(l.y);
			p.associations.push_back(nearest_lm.id_i);
			double e = exp(-(pow((l.x-nearest_lm.x_f)/std_landmark[0],2)+pow((l.y-nearest_lm.y_f)/std_landmark[1],2))/2.0);
			p.weight *= e / (2.0*M_PI*std_landmark[0]*std_landmark[1]);

			//std::cout << "pow((l.x-nearest_lm.x_f)/std_landmark[0],2): "<< pow((l.x-nearest_lm.x_f)/std_landmark[0],2) <<endl;
			//std::cout << "pow((l.y-nearest_lm.y_f)/std_landmark[1],2): "<< pow((l.y-nearest_lm.y_f)/std_landmark[1],2) <<endl;
			//std::cout << "gaussian normalizer: "<<2.0*M_PI*std_landmark[0]*std_landmark[1]<<endl;

		}
		std::cout << "Particle "<<p.id<<" weight: "<<p.weight <<endl;

	}

	std::cout << "Weights updated..." <<endl;

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

		std::vector<double> weights;
		double weight_sum = 0.0;
		// Normalize weights  so total weight is 1.0

		for(Particle p: particles) {
			weight_sum += p.weight;
			weights.push_back(p.weight);
		}
		for(Particle& p : particles) {
			p.weight = p.weight / weight_sum;
		}

	std::default_random_engine generator;
  std::discrete_distribution<> distribution(weights.begin(), weights.end());

	std::vector<Particle> new_particle_set;
	int index;

	for(int i=0; i<particles.size(); i++) {
		index = distribution(generator);
		new_particle_set.push_back(particles[index]);
	}
	particles = new_particle_set;

	std::cout << "Resampling complete..." <<endl;

// Normalize weights  so total weight is 1.0

	weight_sum = 0.0;
	for(Particle p: particles) {
		weight_sum += p.weight;
	}

	for(Particle& p : particles) {
		p.weight = 1.0; //p.weight / weight_sum;
	}

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
