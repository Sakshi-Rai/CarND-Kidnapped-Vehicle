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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cstdlib>
#include <cfloat>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
		num_particles=10;
		default_random_engine gen;
		
	// This line creates a normal (Gaussian) distribution for x,y,theta
		normal_distribution<double> dist_x(x,std[0]);
		normal_distribution<double> dist_y(y,std[1]);
		normal_distribution<double> dist_theta(theta,std[2]);
		
		for(int i=0;i<num_particles;i++){
		Particle p = {};
		p.id=i;
		p.x=dist_x(gen);
		p.y=dist_y(gen);
		p.theta=dist_theta(gen);
		p.weight=1.0;
		particles.push_back(p);
		weights.push_back(p.weight);
		}
		is_initialized=true;
		
		
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	for(int i=0;i<num_particles;i++){
		double in_theta=particles[i].theta;
		if(fabs(yaw_rate)>0.0001){
		particles[i].x += (velocity/yaw_rate)*(sin(in_theta+(yaw_rate*delta_t))-sin(in_theta));
		particles[i].y += (velocity/yaw_rate)*(cos(in_theta)-cos(in_theta+(yaw_rate*delta_t)));
		particles[i].theta += (yaw_rate*delta_t);
		}
		else{
		particles[i].x += (velocity*delta_t*cos(in_theta));
		particles[i].y += (velocity*delta_t*sin(in_theta));
		particles[i].theta = in_theta;
		}
		
		normal_distribution<double> dist_x(0,std_pos[0]);
		normal_distribution<double> dist_y(0,std_pos[1]);
		normal_distribution<double> dist_theta(0,std_pos[2]);
		
		particles[i].x+=dist_x(gen);
		particles[i].y+=dist_y(gen);
		particles[i].theta+=dist_theta(gen);
	}
	//cout << "predicted";
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	//cout << observations.size() << endl;
	std::vector<LandmarkObs> closest_landmark(observations.size());
	double den= (2 * M_PI * std_landmark[0] * std_landmark[1]);
	for(int i=0;i<num_particles;i++)
	{
		double p=1.0;
		double xm ;
		double ym ;
		for(int j=0;j<observations.size();j++)
		{
			
			xm = particles[i].x + (cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y);
			ym = particles[i].y + (sin(particles[i].theta)*observations[j].x) + (cos(particles[i].theta)*observations[j].y);
			//cout << "transformed coordinates" << xm << ym << endl;
			double distance=DBL_MAX;
			for(int k=0;k<map_landmarks.landmark_list.size();k++)
			{
			double value=dist(xm,ym,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
			//cout << "distance" << value << endl;
			if(value<distance)
				{
				//cout << "map landmark" << map_landmarks.landmark_list[k].x_f << map_landmarks.landmark_list[k].y_f << endl;
				closest_landmark[j].x=map_landmarks.landmark_list[k].x_f;
				closest_landmark[j].y=map_landmarks.landmark_list[k].y_f;
				//cout << "closest landmark" << closest_landmark[j].x << "#" <<closest_landmark[j].y << endl;
				distance=value;
				}
			}
			
			double x_term = pow(xm - closest_landmark[j].x, 2) / (2 * pow(std_landmark[0], 2));
			double y_term = pow(ym - closest_landmark[j].y, 2) / (2 * pow(std_landmark[1], 2));
			double w= exp(-(x_term + y_term)) /den;
			//cout << "w" << p <<" "<< w << endl;
			//cout << p << endl;
			p=p*w;
		}
		
		//cout << "final weight" << p << endl;
		particles[i].weight=p;
		weights[i]=p;

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::discrete_distribution<int> ddist(weights.begin(), weights.end());
	default_random_engine gen;
	vector<Particle> resampled_particles;
	for (int i = 0; i < num_particles; i++)
	{
		Particle p={};
		p = particles[ddist(gen)];
		resampled_particles.push_back(p);
	}
	particles = resampled_particles;
	
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
