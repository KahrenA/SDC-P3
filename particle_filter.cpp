/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang;;; + Kahren A 
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

#include "particle_filter.h"

using namespace std;
default_random_engine rand_gen;

//***********************************************************************
void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position 
	// (based on estimates of x, y, theta and their uncertainties from GPS) and 
	// all weights to 1. 
	// Add random Gaussian noise to each particle.

	std::cout << "In ParticleFilter::Init \n";

	// normal dist = (1/(sqrt(std*2PI)) exp (-1/2)( sq((x-ux)/std) )
	std::normal_distribution<double> noise_x (0, std[0]);
	std::normal_distribution<double> noise_y (0, std[1]);
	std::normal_distribution<double> noise_t (0, std[2]);


	num_particles = 10;     // 500 is reasonable # of particles from Sebastian's 

	for (int i=0; i < num_particles; i++)
	{
		Particle particle;
	 	particle.id = i;
		particle.x = x + noise_x(rand_gen);
		particle.y = y + noise_y(rand_gen); 
		particle.theta = theta + noise_t(rand_gen);
		particle.weight = 1.0;
//		particle.associations = null;   // ?
//		particle.sense_x = 0;			// ?
//		particle.sens_y = 0;			// ?

		// Append to particle vector list  
		particles.push_back(particle);	

	}

	is_initialized = true; 
}

//**************************************************************
// 
// prediction
// 
// Update the particle based on reported velocity and yaw_rate
//
void ParticleFilter::prediction(double delta_t, double std_pos[], 
								double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and 
	//std::default_random_engine useful.
	
	std::cout << "In ParticleFilter::prediction\n\n";
	std::cout << "velocity = " << velocity << "\t" << "yaw_rate " << yaw_rate 
					<< "delta_t = " << delta_t << "\n"; 

	double new_x, new_y, new_t; 
	normal_distribution<double> noise_x(0, std_pos[0]);
	normal_distribution<double> noise_y(0, std_pos[1]);
	normal_distribution<double> noise_t(0, std_pos[2]);

	// For all particles 
	for (int i = 0; i < num_particles; i++)
	{
		Particle prtcl = particles[i];

		if(yaw_rate < 0.001)
		{
			new_x = prtcl.x + velocity * delta_t * cos(prtcl.theta);
			new_y = prtcl.y + velocity * delta_t * sin(prtcl.theta); 
			new_t = prtcl.theta;
		}
		else
		{
			new_x = prtcl.x + (velocity/yaw_rate) * 
						 ( sin(prtcl.theta + yaw_rate*delta_t) - sin(prtcl.theta) );
			new_y = prtcl.y + (velocity/yaw_rate) * 
						( cos(prtcl.theta) -  cos(prtcl.theta + yaw_rate*delta_t) );
			new_t = prtcl.theta + yaw_rate * delta_t; 
		}


		prtcl.x = new_x;
		prtcl.y = new_y;
		prtcl.theta = new_t;

		// incorporate noise before saving...
		prtcl.x += noise_x(rand_gen);
		prtcl.y += noise_y(rand_gen);
		prtcl.theta += noise_t(rand_gen);

		std::cout << "new_x = " << prtcl.x << "\t" << "new_y = " << prtcl.y << "\t" 
									<< "new_theta = " << prtcl.theta << "\n";

		particles[i] = prtcl;   // save 
		
	}	// for each particle

	std::cout << "\n\n";
}

//*********************************************************************
//
double distanceBetweenTwoPoints(double x, double y, double a, double b)
{
	return sqrt(pow(x - a, 2) + pow(y - b, 2));
}

//********************************************************************************
//
// dataAssociation 
// 
// Find the predicted measurement that is closest to each observed measurement 
// and assign the observed measurement to this particular landmark.
//
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> inRangeLandmarks, 
									 std::vector<LandmarkObs>& observations) 
{
	
	double min_dist_so_far = std::numeric_limits<double>::max();
	int lm_id=-1;

	for (int i=0; i < observations.size(); i++)
	{
		// find the nearest landmark for this observation
		for(int l=0; l < inRangeLandmarks.size(); l++)
		{
			double distance = 
				distanceBetweenTwoPoints( 
								observations[i].x, observations[i].y,
									inRangeLandmarks[l].x, inRangeLandmarks[l].y);

			if (distance < min_dist_so_far)
			{
				//keep smaller of the two until we get to the smallest
				min_dist_so_far = distance;
				lm_id = inRangeLandmarks[l].id;
			}	
		} // for l

		// Save the Id of landmark closest to observation
		observations[i].id = lm_id;

		cout << "Associated LM ID: " << lm_id << "\t" << "with Observation : " << i << "\n";

	} //for i obs
	
}


//*****************************************************************************
// updateWeights -  
// 
// In this function, we will 
// 1) check and save landmarks that are within sensor range of each particle
// 2) Transform the observations sent by the simulator from vehicle coordinates 
//     to map coordinates 
// 3) Compare and match observations' (measured) distance to landmarks against 
//      particles' (predicted) distance to landmarks. 
// 4) Adjust particle's weight 
//     
//*****************************************************************************
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. 

	std::cout << "In UpdateWeights \n";

	// let's compare particle's distance to landmark with observation's 
	for (int i = 0; i < num_particles; i++)
	{
		Particle prtcl = particles[i];
		cout << "UW : inRangeLandmark id = ";

		// ---------------------------------------------------------------
		// Save landmarks that are within sensor range of the particle 
		// as observations are limited by the sensor range ...
		//vector<single_landmark_s> inRangeLandmarks;
		vector<LandmarkObs> inRangeLandmarks;
		LandmarkObs single_LM;

		for (int l=0; l < map_landmarks.landmark_list.size(); l++)
		{
			if(fabs(map_landmarks.landmark_list[l].x_f - prtcl.x <= sensor_range) && 
				fabs(map_landmarks.landmark_list[l].y_f - prtcl.y <= sensor_range) )
			{
				// Save in inRangeLandmarks
				single_LM.x = map_landmarks.landmark_list[l].x_f;
				single_LM.y = map_landmarks.landmark_list[l].y_f;
				single_LM.id = map_landmarks.landmark_list[l].id_i;
				inRangeLandmarks.push_back( single_LM);

				cout  << single_LM.id << "  ";
			}
		} // for l
		cout << "\n";


		//--------------------------------------------------------------------------
		// First translate observations from vehicle coordinates to map coordinates 
		// Lesson14 sections 13-15 ... note rotation direction is important
		// -------------------------------------------------------------------------
		vector<LandmarkObs> transformed_obs;
		double trans_x, trans_y;
		for (int t =0; t < observations.size(); t++)
		{
			if(prtcl.theta > 0)   // counterclockwise
			{
				trans_x = prtcl.x + cos(prtcl.theta) * observations[t].x
								  - sin(prtcl.theta) * observations[t].y ;
				trans_y = prtcl.y + sin(prtcl.theta) * observations[t].x
								  + cos(prtcl.theta) * observations[t].y; 
			}
			else			// Clockwise
			{
				trans_x = prtcl.x + cos(prtcl.theta) * observations[t].x
								  + sin(prtcl.theta) * observations[t].y ;
				trans_y = prtcl.y - sin(prtcl.theta) * observations[t].x
								  + cos(prtcl.theta) * observations[t].y; 
			}
			LandmarkObs single_translated_observation;
			single_translated_observation.x = trans_x;
			single_translated_observation.y = trans_y;
			single_translated_observation.id = observations[t].id;

			cout << "OBSs: trans_x = " << trans_x << "\t" << "trans_y = " << trans_y << "\n";
			transformed_obs.push_back(single_translated_observation);
		}

		//---------------------------------------------------------------------
		// Associate measured (and translated) observations with inRangelandmarks 
		// by updating the id of the observation with the landmark id
		//---------------------------------------------------------------------
		dataAssociation( inRangeLandmarks, transformed_obs);



		//-------------------------------------------------------------------------
		// Update weight of this particle using multi-variate_normal_dist Lesson14
		// sections 14-18
		// # calculate normalization term
		// gauss_norm= (1/(2 * np.pi * sig_x * sig_y))
		//
		// # calculate exponent
		// exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + 
		//					((y_obs - mu_y)**2)/(2 * sig_y**2)
		//
		// # calculate weight using normalization terms and exponent
		// weight= gauss_norm * math.exp(-exponent)
		//-------------------------------------------------------------------------
		particles[i].weight = 1.0;

		// From the set of observations which now have landmarkIDs, 
		for(int t=0; t < transformed_obs.size(); t++)
		{
			for (int l=0; l < inRangeLandmarks.size(); l++)
			{
				//let find the landmark that is seen by this observation
				if(transformed_obs[t].id == inRangeLandmarks[l].id)
				{
					// Update weight of particle based on formula in comment above
					double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
					// cout << "gauss_norm = " << gauss_norm << "\t";

					double xexp = pow( (transformed_obs[t].x) - inRangeLandmarks[l].x, 2) / 
														(2*pow(std_landmark[0], 2));
										
					cout << "inRangeLandmarks[l].x = " << inRangeLandmarks[l].x  << "\t" 
						 << "transformed_obs[t].x  = " << transformed_obs[t].x << "\t" 
						 << "xexp = " << xexp << "\n";

					double yexp = pow( (transformed_obs[t].y - inRangeLandmarks[l].y), 2) / 
														(2*pow(std_landmark[1], 2));
					cout << "inRangeLandmarks[l].y = " << inRangeLandmarks[l].y << "\t" 
						 << "transformed_obs[t].y  = " << transformed_obs[t].y << "\t" 
						 << "yexp = " << yexp << "\n";


					double weight = gauss_norm * exp(-(xexp +yexp));

					cout << "weight = " << weight << "\t" << "particle[i].weight = " << particles[i].weight << "\n";

					particles[i].weight *= weight; 

					cout << "New weight = " << particles[i].weight << "\n";

					break;
				}
			} // for l

		} // for t 

 	} // i num_particles

}

//******************************************************************
// 
// resample
// 
// Based on Sebastian's wheel 
// 
// index = random # between 0 -> N-1
// Beta = 0
// Beta = Beta + U (continuous space from 0 --> 2 * max weight)
// if w[index] < Beta
// 	    Beta = Beta - weight[index]
// 		index = index + 1
// else
//		pick particle at index
//
//******************************************************************	
void ParticleFilter::resample() 
{
	
	std::cout << "In resample \n";
	
	double beta = 0.0;
	vector <double> weights;
	
	// determine max weight
	for (int i=0; i < num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}
	double max_w = *max_element(weights.begin(), weights.end()); 
	cout << "Max weight: "<<max_w << endl;

	//--------------------------
	// Implement the wheel algo
	//--------------------------
	std::default_random_engine eng;
	std::uniform_int_distribution<int> dist{0, num_particles-1};
    unsigned int index = dist(eng);
	vector<Particle> resampled_particles;	

	for (int i=0; i < num_particles; i++)
	{
		beta = beta + 2 * max_w;
		while (weights[index] < beta) 
		{
			beta = beta - weights[index];
			index = ( index + 1 ) % num_particles;
		}
		
		// take the particle into new set 
		resampled_particles.push_back(particles[index]);
	} 
	
	// Replace old set of particles with the new set 
	particles = resampled_particles;
}

//**********************************************************************************************
Particle ParticleFilter::SetAssociations(Particle particle, 
											std::vector<int> associations, 
											std::vector<double> sense_x, 
											std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) 
	// world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

//**********************************************************
string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));

    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

//**********************************************************
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));

    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

//**********************************************************
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));

    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
