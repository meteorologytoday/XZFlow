#define DEBUG_MODE

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#define _USE_MATH_DEFINES
#include <cmath>



#ifdef DEBUG_MODE
	#define __DBG__ 
#else
	#define __DBG__ //
#endif

#define JACOBIAN_FAILED 11



using namespace std;

typedef float T;

const T kappa = 2.0 / 7.0;
const T cp_over_cv = 7.0 / 3.0;
const T R_ps_kappa = 287.0 / pow(101300, 2.0 / 7.0);
const T g0 = 9.8;


const int xpts = 256;
const int zpts = 256;
const int pts = xpts * zpts;

const T Lx = 10000.0f;
const T Lz = 10000.0f;
const T nu = 0.1;
const T theta_mean = 300.0f;


T dx = Lx / (xpts - 1);
T dz = Lz / (zpts - 1);

T *zeta, *theta, *rho, *psi, *u, *w, *gdx_zeta, *gdz_zeta, \
	*gdx_theta, *gdz_theta, *lap_zeta, *torq,          \
	*zeta_k1,  *zeta_k2,  *zeta_k3,  *zeta_k4,         \
	*theta_k1, *theta_k2, *theta_k3, *theta_k4,        \
	*tmp, *zeta_tmp, *rho_tmp, *theta_tmp;

const int total_steps = 3000;
const int record_step = 100;
const T dt = 3;



const int max_iter = 5000;
const T max_err = 0.01 * dx;

inline int I(int i, int j) {
	return j * xpts + i;
}

void init();
void gradx(T *, T *);
void gradz(T *, T *);
void evolve(T *, T *, T, T*);
int jacobi_relaxation(T * input, T * output, T * tmp, T max_err, int max_iter);
void laplacian(T * input, T * output);
int output_data(T *, const char *);

void init() {
	// allocate arrays
	zeta       = new T [pts];
	theta      = new T [pts];
	rho        = new T [pts];
	psi        = new T [pts];
	u          = new T [pts];
	w          = new T [pts];
	gdx_zeta   = new T [pts];
	gdz_zeta   = new T [pts];
	gdx_theta  = new T [pts];
	gdz_theta  = new T [pts];
	lap_zeta   = new T [pts];
	torq       = new T [pts];
	tmp        = new T [pts];
	zeta_tmp   = new T [pts];
	theta_tmp  = new T [pts];

	// used for RK4	
	zeta_k1    = new T [pts];
	zeta_k2    = new T [pts];
	zeta_k3    = new T [pts];
	zeta_k4    = new T [pts];
	
	theta_k1   = new T [pts];
	theta_k2   = new T [pts];
	theta_k3   = new T [pts];
	theta_k4   = new T [pts];

	// init field
	float x, z, r2;
	float sigma = 1000.0;
	float constant = 3000.0 / sigma / pow(2.0 * 3.14, 0.5);
	float cx = 5000.0, cz = 3000.0;
	for(int j=0; j < zpts; ++j) {
		z = dz * j;
		for(int i=0; i < xpts; ++i) {
			x = dx * i;
			r2 = pow(x - cx, 2.0) + pow(z - cz, 2.0);
			
			theta[I(i,j)] = constant * exp(- r2 / pow(sigma, 2.0)); 
		}
	}
}



void cal_change(T * zeta, T * theta, T * zeta_change, T * theta_change) {

	// calculate torq
	gradx(theta, torq);

	for(int i = 0; i < pts; ++i) {
		torq[i] *= g0 / theta_mean;
	}
	output_data(torq, "torq.bin");
	
	int flag = 0;
	// invert wind from zeta
	if((flag = jacobi_relaxation(zeta, psi, tmp, max_err, max_iter)) != 0) {
		fprintf(stderr, "[Error code %d]: Relaxation cannot converge within [%d] iterations.\n", flag, max_iter);
		throw JACOBIAN_FAILED;
	}

	// calculate u, w
	gradz(psi, u); for(int i =0; i < pts; ++i) { u[i] *= -1.0f; } 
	gradx(psi, w);

	output_data(u, "u.bin");
	output_data(w, "w.bin");
	// calculate grad vorticity
	gradx(zeta, gdx_zeta);
	gradz(zeta, gdz_zeta);

	// calculate lap vort
	laplacian(zeta, lap_zeta);

	for(int i = 0; i < pts; ++i) {
		zeta_change[i] = - u[i] * gdx_zeta[i] - w[i] * gdz_zeta[i] + nu * lap_zeta[i] + torq[i];
	}
	
	// calculate theta part
	gradx(theta, gdx_theta);
	gradz(theta, gdz_theta);
	for(int i = 0; i < pts; ++i) {
		theta_change[i] = - u[i] * gdx_theta[i] - w[i] * gdz_theta[i];
	}
	output_data(theta_change, "theta_change.bin");	
}


int main() {

	init();
	int record_flag;
	char filename[256];
	for(int step = 0; step < total_steps; ++ step) {
		printf("# Step %d, time = %.2f", step, step * dt);
		if( (record_flag = ((step % record_step) == 0)) ) { printf(", record now!");}
		printf("\n");

		if(record_flag) {
			sprintf(filename, "output/zeta_step_%d.bin", step);
			output_data(zeta, filename);
			
			sprintf(filename, "output/theta_step_%d.bin", step);
			output_data(theta, filename);
		}

		for(int rk = 0 ; rk < 4; ++rk) {

			switch (rk) {
				case 0:
					cal_change(zeta, theta, zeta_k1, theta_k1);
					break;
				case 1:
					evolve(zeta,  zeta_k1,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k1, dt/2.0, theta_tmp);
					
					cal_change(zeta_tmp, theta_tmp, zeta_k2, theta_k2);

					break;
				case 2:
					evolve(zeta,  zeta_k2,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k2, dt/2.0, theta_tmp);
					
					cal_change(zeta_tmp, theta_tmp, zeta_k3, theta_k3);
					break;
				case 3:
					evolve(zeta,  zeta_k3,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k3, dt/2.0, theta_tmp);
					
					cal_change(zeta_tmp, theta_tmp, zeta_k4, theta_k4);
					break;
			}
		}

		for(int i=0; i<pts; ++i) {
			zeta[i] += (zeta_k1[i] + 2.0 * zeta_k2[i] + 2.0 * zeta_k3[i] + zeta_k4[i]) / 6.0;
			theta[i] += (theta_k1[i] + 2.0 * theta_k2[i] + 2.0 * theta_k3[i] + theta_k4[i]) / 6.0;
		}
	}
}

const T d2x = 2 * dx;
void gradx(T * input, T * output) {
	for(int j=0; j<zpts; ++j) {
		for(int i=0; i<xpts; ++i) {
			if(i == 0) {
				output[I(i,j)] = (input[I(i+1,j)] - input[I(i,j)]) / dx;
			} else if(i == xpts-1) {
				output[I(i,j)] = (input[I(i,j)] - input[I(i-1,j)]) / dx;
			} else {
				output[I(i,j)] = (input[I(i+1, j)] - input[I(i-1,j)]) / d2x;
			}
		}
	}
}

const T d2z = 2 * dz;
void gradz(T * input, T * output) {
	for(int j=0; j<zpts; ++j) {
		for(int i=0; i<xpts; ++i) {
			if(j == 0) {
				output[I(i,j)] = (input[I(i,j+1)] - input[I(i,j)]) / dz;
			} else if(j == zpts-1) {
				output[I(i,j)] = (input[I(i,j)] - input[I(i,j-1)]) / dz;
			} else {
				output[I(i,j)] = (input[I(i, j+1)] - input[I(i,j-1)]) / d2z;
			}
		}
	}
}

void evolve(T * f, T * df, T dt, T * o_f) {
	for(int i=0; i < pts; ++i){
		o_f[i] = f[i] + df[i] * dt;
	}
}

int jacobi_relaxation(T * input, T * output, T * tmp, T max_err, int max_iter) {
	/*
	 * Initial guess is stored in output
	 * Boundary condition is specified on the margin of input
	 */
	T err = 0, r;
	T gamma = 1 / (2 * ( 1/(dx*dx) + 1/(dz*dz) ) );
	int avg_pts = (xpts-2) * (zpts-2);
	for(int iter = 0; iter < max_iter; ++iter) {
		
		err = 0;
		// calculate nth step zeta field
		laplacian(output, tmp);

		// tune
		for(int i = 1; i < xpts-1; ++i) {
			for(int j = 1; j < zpts-1; ++j) {
				r = tmp[I(i,j)] - input[I(i,j)];
				err += r * r;
				output[I(i,j)] += gamma * r;
			}
		}

		err /= avg_pts;
		fprintf(stderr,"[%d] Err avg: %e\n",iter, err);

		if(err < max_err)
			return 0;
	}
	
	return -1;
}

void laplacian(T * input, T * output) {
	T dxdx = dx * dx;
	T dzdz = dz * dz;
	int centx, centz;
	for(int i = 0; i < xpts; ++i) {
		for(int j = 0; j < zpts; ++j) {
			centx = i; centz = j;
			if(i == 0)
				++centx;
			else if(i == xpts-1)
				--centx;

			if(j == 0)
				++centz;
			else if(j == zpts-1)
				--centz;

			output[I(i,j)] = (input[I(centx+1,centz)] + input[I(centx-1,centz)] - 2*input[I(centx,centz)]) / dxdx \
						  +  (input[I(centx,centz+1)] + input[I(centx,centz-1)] - 2*input[I(centx,centz)]) / dzdz;
		}
	}

}

int output_data(T * dat, const char * filename) {
	FILE * file = fopen(filename,"wb");
	fwrite(dat, sizeof(T), pts, file);
	fclose(file);
	return 0;
}
