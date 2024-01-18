#include "simulation.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

double randn(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  call = !call;
  return (mu + sigma * (double) X1);
}

Particle* initialize_particle(double mean_mass,double var_mass,double mean_pos,double var_pos){
    Particle* p = malloc(sizeof(Particle));
    p->m = randn(mean_mass,var_mass);
    for(int i = 0;i<dimensions;i++){
        p->r[i] = randn(mean_pos,var_pos);
        p->u[i] = p->v[i] = 0;
    }
    return p;
}

Particle** initialize_particle_system(double mean_mass,double var_mass,double mean_pos,double var_pos,int n){
    Particle **p = malloc(sizeof(Particle)*n);
    for(int i = 0;i<n;i++){
        p[i] = initialize_particle(mean_mass,var_mass,mean_pos,var_pos);
    }
    return p;
}

void collision(Particle *p1,Particle *p2){
    if(&p1 == NULL || &p2 == NULL) return;
    Particle p;
    p.m = p1->m + p2->m;
    for(int i = 0;i<dimensions;i++){
        p.r[i] = p1->m*p1->r[i] + p2->m*p2->r[i];
        p.u[i] = p1->m*p1->u[i] + p2->m*p2->u[i];
        p.v[i] = p1->m*p1->v[i] + p2->m*p2->v[i];
    }
    copy_params(&p,p1);
    free(p2);
    p2 = NULL;
}


Particle* update_particle_params(Particle **p,int particle_num){
    if(p[particle_num] == NULL) return;
    int n = sizeof(p)/sizeof(Particle);
    double g[dimensions] = {0.0};
    for(int i = 0;n;i++){
        if(i == particle_num) continue;
        if(p[i] == NULL) continue;
        double r = calc_euclidean_dist_cubed(p[i],p[particle_num]);
        for(int j = 0;j<dimensions;j++){
            g[j] = (G*(p[i]->m)*(p[i]->r[j] - p[particle_num]->r[j]))/r;
        }
    }
    Particle *temp_p;
    for(int i = 0;i<dimensions;i++){
        temp_p->u[i] = p[particle_num]->u[i];
        temp_p->v[i] = p[particle_num]->v[i];
        temp_p->r[i] = p[particle_num]->r[i];
    }
    return temp_p;
}

void update_system(Particle **p){
    int n = sizeof(p)/sizeof(Particle);
    Particle** p_temp = initialize_particle_system(0,0,0,0,n);
    for(int i = 0;i<n;i++){
        p_temp[i] = update_particle_params(p,i);
    }
    copy_particle_system(p_temp,p);
}


double calc_euclidean_dist_cubed(Particle *p1,Particle *p2){
    if(p1 == NULL || p2 == NULL) return;
    double ans = 0.0;
    for(int i = 0;i<dimensions;i++){
        ans += pow((p1->r[i] - p2->r[i]),2);
    }
    return pow(sqrt(ans),2);
}

void copy_params(Particle *p1,Particle *p2){
    if(p1 == NULL || p2 == NULL) return;
    p2->m = p1->m;
    for(int i = 0;i<dimensions;i++){
        p2->r[i] = p1->r[i];
        p2->u[i] = p1->u[i];
        p2->v[i] = p1->v[i];
    }
    
}

void copy_particle_system(Particle **p1,Particle **p2){
    int n = sizeof(p1)/sizeof(Particle);
    for(int i = 0;i<n;i++){
        copy_params(p1,p2);
    }
}

void delete_particle_system(Particle **p){
    if(p == NULL) return;
    int n = sizeof(p)/sizeof(Particle);
    for(int i = 0;i<n;i++){
        free(p[i]);
        p[i] = NULL;
    }
    free(p);
    p = NULL;
}

void print_particle_info(Particle *p){
    if(p == NULL) return;
    printf("Mass of the particle is: %f\n",p->m);
    printf("Position of the particle is: [");
    for(int i = 0;i<dimensions;i++){
        printf("%f ",p->r[i]);
    }
    printf("]\n");
    printf("Initial Velocity of the particle is: [");
    for(int i = 0;i<dimensions;i++){
        printf("%f ",p->r[i]);
    }
    printf("]\n");
    printf("Initial Velocity of the particle is: [");
    for(int i = 0;i<dimensions;i++){
        printf("%f ",p->u[i]);
    }
    printf("]\n");
    printf("Final Velocity of the particle is: [");
    for(int i = 0;i<dimensions;i++){
        printf("%f ",p->v[i]);
    }
    printf("]\n");
}

int main(){
    printf("%ld\n",1<<30);
    printf("%ld\n",((1ll)*sizeof(Particle)*(1ll)*INT_MAX)/(1<<30));
    // int n = INT_MAX;
    // printf("Dimensions: %d\n",dimensions);
    // Particle **p = initialize_particle_system(1,1,0,1,n);
    // for(int i = 0;i<n;i+=INT_MAX/100) print_particle_info(p[i]);
    return 0;
}
