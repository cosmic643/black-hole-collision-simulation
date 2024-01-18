#pragma once

#define dimensions 3
const double merger_radius = 0.0;
const double iteration_time_scale = 0.0;
const double G = 6.6743/1e11;
typedef struct{
    double m;
    double r[dimensions];
    double u[dimensions];
    double v[dimensions];
}Particle;

Particle* initialize_particle(double mean_mass,double var_mass,double mean_pos,double var_pos);
Particle** initialize_particle_system(double mean_mass,double var_mass,double mean_pos,double var_pos,int n);
Particle* update_particle_params(Particle **p,int particle_num);
void delete_particle_system(Particle **p); //utility function
void copy_params(Particle *p1,Particle *p2); //copies p1's data in p2
double calc_euclidean_dist_cubed(Particle *p1,Particle *p2); //will calculate the 3/2th power of the eucledian distance
void print_particle_info(Particle *p); //utility functions
double randn(double mu, double sigma);
void copy_particle_system(Particle **p1,Particle **p2);