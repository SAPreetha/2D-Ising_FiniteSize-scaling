//
//  Preetha_Project2.cpp
//  
//
//  Created by Preetha Saha on 4/6/17.
//
//


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>


using namespace std;

const int N = 4;// linear lattice size

const int Ns = N * N;   // number of spins
const int N_nn = 4; // No of Nearest neighbour of each lattice point

const double J = 1; //ferromagnetic interaction in nearest neighbours
float T;

const int N_eq=5000; //Thermalization steps
const int N_MC=100000; // No of monte carlo steps
int acc_count;//number of spins flipped



double pow(double x, double y);
double fabs(double x);
//double ln(double x);
double log(double x);



class lattice{ //Defining Lattice class to keep spin variables
public:
    int idx; //index of lattice point in stored 1D array
    int x; //x position of lattiice point
    int y; //y position of lattice point
    int Sz; // Sz spin value, can be +1 or -1
    
    
    lattice *nn[N_nn];
    
};

lattice ising_spin[Ns]; // lattice objet to store spin values

//return mod value to take care of pereodic boundary condition
int mod(int x, int m) {
    if (x>=0 && x<m)
        return x;
    else if (x<0)
        return m-1-mod(-1-x,m);
    else
        return x%m;
}

// 2D variable stored in 1D array, indexed accordingly
int index_site(int x, int y){
    int xm=mod(x,N);
    int ym=mod(y,N);
    
    return (xm*N)+ym;
}

//lattice position variables defined

void set_lattice_variables(void){
    for(int i=0;i<Ns;i++){
        for(int j=0;j<Ns;j++){
            int idx=index_site(i,j);
            ising_spin[idx].idx=idx;
            ising_spin[idx].x=i;
            ising_spin[idx].y=j;
            
        }
    }
}
//nearest neighbour indices
void set_nn(void){
    for(int i=0;i<Ns;i++){
        int x=ising_spin[i].x;
        int y=ising_spin[i].y;
        int j;
        
        j=index_site(x+1,y);
        ising_spin[i].nn[0]=&ising_spin[j];
        
        j=index_site(x-1,y);
        ising_spin[i].nn[1]=&ising_spin[j];
        
        j=index_site(x,y+1);
        ising_spin[i].nn[2]=&ising_spin[j];
        
        j=index_site(x,y-1);
        ising_spin[i].nn[3]=&ising_spin[j];
        
        
    }
}
//intial lattice spins randomly initiated
void initiate_lattice(void){
    for(int i=0; i<Ns; i++){
        float temp=(float)rand()/RAND_MAX;
        ising_spin[i].Sz = (temp > 0.5) ? +1 : -1;
        printf("%d\n",ising_spin[i].Sz);
        
    }
    
}
// Change in Energy of lattice on flipping a spin

double delta_Energy(int random_pos){
    
    int sum=0;
    
    for(int n=0;n<N_nn;n++){
        
        sum=sum+ising_spin[random_pos].nn[n]->Sz;
    }
    
    
    double dE=-J * (-2*ising_spin[random_pos].Sz) * sum;
    
    return dE;
}

//function of update the spinflip according to dE conditions, which returns the no of spins flipped if 0 or 1

double spinflip_update(float T){
    int random_pos=rand()%Ns;
    double delta_E = delta_Energy(random_pos);
    
    double acc_count=0;
    
    if (delta_E <0){
        ising_spin[random_pos].Sz*=-1;
        acc_count += 1;
        
        
    }
    if (delta_E >= 0){
        float rand_n=(float)rand()/RAND_MAX;
        double prob=exp(-delta_E/T);
        if(prob>rand_n){
            ising_spin[random_pos].Sz*=-1;
            acc_count += 1;
            
            
            
        }
    }
    return acc_count;
    
}

// Calculating Energy of the lattice
double Energy(void){
    double Ene=0;
    for(int i=0;i<Ns;i++){
        int sum=0;
        for(int n=0;n<N_nn;n++){
            sum += ising_spin[i].nn[n]->Sz;
        }
        Ene += -J * ising_spin[i].Sz * sum;
        
        
    }
    return (double)Ene/2;
}

// Calculating Magnetization of the lattice
int Magnetisation(void){
    int Mag=0;
    for(int i=0;i<Ns;i++){
        Mag = Mag + ising_spin[i].Sz;
        
    }
    
    return Mag;
    
}
//Metropolis algorithm= Attempting Ns spin flips and returning new lattice and total no. of spins flipped.

double metropolis(float T){
    double C=0;
    double Acc_rate;
    for(int i=0;i<Ns;i++){
        double Ene=0;
        C += spinflip_update(T);
        
    }
    
    
    Acc_rate=(double)C/Ns;
    return Acc_rate;//Acceptance of spin flip rate over every Monte carlo step
    
}
// Getting the system in Equilibrium before taking measurements

void Thermal_eq_sweep(void){
    for(int i=0;i<N_eq;i++){
        double R;
        
        R= metropolis(T);
    }
}

//Looping over Teperature
void Measurement(){
    
    FILE *g = fopen("datN.txt", "w+");
    
    for(double T=.4;T<4;T+=.02){
        
        
        double Enesum=0;
        double Magsum=0;
        double Mabs_sum=0;
        
        double Mag_sq_sum=0;

        double Ene_sq_sum=0;
        double count=0;
        double M4sum=0;
        double Rsum=0;
        

        
        for(int i=0;i<N_MC;i++){
            double R=0;
            R=metropolis(T);// One Montecarlo step
            if(i%100==0){
                count +=1;
                Rsum += R;
                
                double Ene =  Energy();
                double Mag = Magnetisation();
                double M_abs=fabs(Mag);
                
                
                Enesum = Enesum+Ene;
                Magsum = Magsum+Mag;
                Mabs_sum += M_abs;
                Ene_sq_sum += Ene*Ene;
                Mag_sq_sum += Mag*Mag;
                M4sum += pow(Mag,4.0);

                
                
                
            }
        }
        double Ravg=Rsum/count;
        double E_avg=Enesum/(count*Ns);
        double m=Mabs_sum/(count*Ns);
        double m_scaled= m*(pow((double)N,(1.0/8)));//beta=1/8
        
        double B4=1-(M4sum/count)/(3*pow((Mag_sq_sum/count),2.0));
        
        double Chi=((Mag_sq_sum/count)-((Mabs_sum/count)*(Mabs_sum/count)))/(T*Ns);
        double Chi_scaled=Chi*(pow((double)N,(-7.0/4)));//gamma 7/8
        
        double Cv=((Ene_sq_sum/(count)-((Enesum/count)*(Enesum/count)))/(T*T))/Ns;
        double Cv_scaled=Cv/log(N);
        double t_scaled=(double)(T-2.27)*N;

        //printf("Final Enesum is %lf\n",Enesum/(count*Ns));
        fprintf(g,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",T,t_scaled,B4,m,m_scaled,Chi,Chi_scaled,Cv,Cv_scaled);
        
    }
}


int main(){
    
    srand(time(NULL));
    set_lattice_variables();
    set_nn();
    initiate_lattice();
    Thermal_eq_sweep();
    Measurement();
    
    return 0;
    
}























