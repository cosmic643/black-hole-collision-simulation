#include"simulation.h"
#include"pbPlots.h"
#include"supportLib.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

double TIME_SCALE = 0.01;

double randn(double mu, double sigma)
{
    double u1 = rand() / (double)RAND_MAX;
    double u2 = rand() / (double)RAND_MAX;

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

    // Apply mean and standard deviation
    return mu + sigma * z0;
}

Particle* initialize_particle(double mean_mass,double var_mass,double mean_pos,double var_pos){
    Particle* p = malloc(sizeof(Particle));
    p->m = randn(mean_mass,var_mass);
    for(int i = 0;i<dimensions;i++){
        p->r[i] = randn(mean_pos,var_pos);
        p->u[i] = 0;
        p->v[i] = 0;
    }
    return p;
}

Particle** initialize_particle_system(double mean_mass,double var_mass,double mean_pos,double var_pos,int n){
    Particle **p = malloc(sizeof(Particle*)*n);
    for(int i = 0;i<n;i++){
        p[i] = initialize_particle(mean_mass,var_mass,mean_pos,var_pos);
    }
    return p;
}

void collision(Particle **p,int particle_num,int n){
    if(p == NULL || p[particle_num] == NULL) return;
    int arr[n];
    int k = 0;
    for(int i = 0;i<n;i++){
        if(i == particle_num) continue;
        if(p[i] == NULL) continue;
        double r = calc_euclidean_dist(p[particle_num],p[i]);
        if(r<=merger_radius) arr[k++] = i;
    }
    if(k == 0) return;
    printf("\nNumber of Collisions: %d\n",k);
    double net_mass = p[particle_num]->m,r[dimensions],u[dimensions];
    for(int i = 0;i<dimensions;i++){
        r[i] = 0.0;
        u[i] = 0.0;
    }
    for(int i = 0;i<k;i++){
        if(i == particle_num) continue; 
        net_mass += p[arr[i]]->m;
        for(int j = 0;j<dimensions;j++){
            r[j] += p[arr[i]]->r[j]*p[arr[i]]->m;
            u[j] += p[arr[i]]->u[j]*p[arr[i]]->m;
        }
        p[arr[i]] = NULL;
    }
    p[particle_num]->m = net_mass;
    for(int i = 0;i<dimensions;i++){
        p[particle_num]->r[i] = r[i]/net_mass;
        p[particle_num]->u[i] = u[i]/net_mass;
        p[particle_num]->v[i] = 0;
    }
}

Particle* update_particle_params(Particle **p,int particle_num,int n,int t){
    if(p[particle_num] == NULL) return NULL;
    double g[dimensions];
    for(int i = 0;i<dimensions;i++) g[i] = 0.0;
    for(int i = 0;i<n;i++){
        if(i == particle_num) continue;
        if(p[i] == NULL) continue;
        double r = calc_euclidean_dist(p[i],p[particle_num]);
        for(int j = 0;j<dimensions;j++){
            g[j] += ((p[i]->m)*(p[i]->r[j] - p[particle_num]->r[j]))/(r*r*r);
        }
    }
    Particle *temp_p = initialize_particle(0,0,0,0);
    for(int i = 0;i<dimensions;i++){
        temp_p->u[i] = p[particle_num]->v[i];
        temp_p->r[i] = p[particle_num]->r[i] + (p[particle_num]->u[i] * TIME_SCALE) + (0.5*g[i]*TIME_SCALE*TIME_SCALE);
        temp_p->v[i] = p[particle_num]->u[i] + g[i]*TIME_SCALE;
    }
    temp_p->m = p[particle_num]->m;
    return temp_p;
}

void update_system(Particle **p,int n,int t){
    Particle **p_temp = initialize_particle_system(0,0,0,0,n);
    for(int i = 0;i<n;i++){
        p_temp[i] = update_particle_params(p,i,n,t);
    }
    copy_particle_system(p_temp,p,n);
    delete_particle_system(p_temp,n);
    for(int i = 0;i<n;i++){
        if(p[i] == NULL) continue;
        //printf("%d\n",i);
        collision(p,i,n);
    }
}


double calc_euclidean_dist(Particle *p1,Particle *p2){
    if(p1 == NULL || p2 == NULL) return 0.0; //wrong
    double ans = 0.0;
    for(int i = 0;i<dimensions;i++){
        ans += pow((p1->r[i] - p2->r[i]),2);
    }
    return sqrt(ans);
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

void copy_particle_system(Particle **p1,Particle **p2,int n){
    for(int i = 0;i<n;i++){
        copy_params(p1[i],p2[i]);
    }
}

void delete_particle_system(Particle **p,int n){
    if(p == NULL) return;
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
        printf("%f ",p->u[i]);
    }
    printf("]\n");
    printf("Final Velocity of the particle is: [");
    for(int i = 0;i<dimensions;i++){
        printf("%f ",p->v[i]);
    }
    printf("]\n\n");
}

void print_system_info(Particle **p,int n){
    printf("\nNumber of Particles in the System: %d\n",n);
    for(int i = 0;i<n;i++) {
        if(p[i] == NULL)continue;
        print_particle_info(p[i]);
    }
}

int count_particles(Particle **p,int n){
    int particle_count = 0;
    for(int i = 0;i<n;i++){
        if(p[i] != NULL){
            particle_count++;
        }
        
    }

    return particle_count;
}

int main(){
    int n = 10;
    int itr = 1000;
    // printf("\t\t\t\t\t\tBlack Hole Collosion Simulation\n");
    // printf("Enter the number of Black holes in Simulation: ");
    // scanf("%d",&n);
    // printf("Enter the number of iterations to be done: ");
    // scanf("%d",&itr);
    // printf("Enter the time to be elapsed in iteration: ");
    // scanf("%f",&TIME_SCALE);
    // printf("Number of Black Holes: %d\n",n);
    // printf("Number of iterations: %d\n",itr);
    // printf("Time elapsed in one iteration: %f\n",TIME_SCALE);
    //printf("\nDimensions: %d\n",dimensions);

    double time[itr],pno[itr];
    Particle **p = initialize_particle_system(10,2,0,1,n);
    print_system_info(p,n);
    printf("----------------\n");
    for(int i = 0;i<itr;i++){
        int temp = count_particles(p,n);
        update_system(p,n,i);
        time[i] = i;
        pno[i] = temp;
        //print_system_info(p,n);
    }
    // for(int i = 0;i<itr;i++){
    //     printf("Time: %f  Number of Particles: %f\n",time[i],pno[i]);
    // }

    //code for storing data in a text file used for graphing purpose
    FILE *file1 = fopen("time.txt", "w");
    for (int i = 0; i < itr; i++) {
        fprintf(file1, "%.2lf\n", time[i]);
    }
    fclose(file1);
    FILE *file2 = fopen("pno.txt", "w");
    for (int i = 0; i < itr; i++) {
        fprintf(file2, "%.2lf\n", pno[i]);
    }
    fclose(file2);
    // RGBABitmapImageReference *imageRef = CreateRGBABitmapImageReference();
    // DrawScatterPlot(imageRef,600,400,time,2,pno,2);
    // size_t length;
    // double *pngData = ConvertToPNG(&length,imageRef->image);
    // WriteToFile(pngData,length,"plot.png");
    return 0;
}