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

const int xpts = 256;
const int zpts = 256;
const int pts = xpts * zpts;

const T Lx = 600000.0f;
const T Lz = 600000.0f;
const T nu = 0.1;
T dx = Lx / (xpts - 1);
T dz = Lz / (zpts - 1);

T *zeta, *theta, *rho, *psi, *u, *w, *gdx_zeta, *gdz_zeta, *gdx_rho, *gdz_rho, \
	*gdx_theta, *gdz_theta, *lap_zeta, *PGF_x, *PGF_z, *torq_1, *torq_2, *total_torq, *pre, \
	*zeta_k1,  *zeta_k2,  *zeta_k3,  *zeta_k4, \
	*theta_k1, *theta_k2, *theta_k3, *theta_k4, \
	*rho_k1,   *rho_k2,   *rho_k3,   *rho_k4, \
	
	*tmp, *zeta_tmp, *rho_tmp, *theta_tmp;

const int total_steps = 1000;
const int record_step = 10;
const T dt = 3;

const int max_iter = 5000;
const T max_err = 0.001;

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
	pre        = new T [pts];
	u          = new T [pts];
	w          = new T [pts];
	gdx_zeta   = new T [pts];
	gdz_zeta   = new T [pts];
	gdx_rho    = new T [pts];
	gdz_rho    = new T [pts];
	gdx_theta  = new T [pts];
	gdz_theta  = new T [pts];
	lap_zeta   = new T [pts];
	PGF_x      = new T [pts];
	PGF_z      = new T [pts];
	torq_1     = new T [pts];
	torq_2     = new T [pts];
	total_torq = new T [pts];
	tmp        = new T [pts];
	zeta_tmp   = new T [pts];
	theta_tmp  = new T [pts];
	rho_tmp    = new T [pts];

	// used for RK4	
	zeta_k1    = new T [pts];
	zeta_k2    = new T [pts];
	zeta_k3    = new T [pts];
	zeta_k4    = new T [pts];
	
	theta_k1   = new T [pts];
	theta_k2   = new T [pts];
	theta_k3   = new T [pts];
	theta_k4   = new T [pts];

	rho_k1     = new T [pts];
	rho_k2     = new T [pts];
	rho_k3     = new T [pts];
	rho_k4     = new T [pts];

	// assign profile
	for(int j = 0; j < zpts; ++j) {
		for(int i=0; i < xpts; ++i) {
			theta[I(i,j)] = 300.0f;
			pre[I(i,j)] = 
		}
	}
}



void cal_change(T * zeta, T * theta, T * rho, T * zeta_change, T * theta_change, T * rho_change) {
	// calculate pressure
	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			pre[I(i,j)] = pow(rho[I(i,j)] * theta[I(i,j)] * R_ps_kappa, cp_over_cv);
		}
	}

	// calculate PGF
	gradx(pre, PGF_x);
	gradz(pre, PGF_z);

	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			PGF_x[I(i,j)] /= - rho[I(i,j)];
			PGF_z[I(i,j)] /= - rho[I(i,j)];
		}
	}

	// calculate torq
	gradz(PGF_x, torq_1);
	gradx(PGF_z, torq_2);
	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			total_torq[I(i,j)] = torq_1[I(i,j)] - torq_2[I(i,j)];
		}
	}

	int flag = 0;
	// invert wind from zeta
	if((flag = jacobi_relaxation(zeta, psi, tmp, max_err, max_iter)) != 0) {
		fprintf(stderr, "[Error code %d]: Relaxation cannot converge within [%d] iterations.\n", flag, max_iter);
		throw JACOBIAN_FAILED;
	}

	// calculate u, w
	gradz(psi, u); for(int i =0; i < pts; ++i) { u[i] *= -1.0f; } 
	gradx(psi, w);

	// calculate grad vorticity
	gradx(zeta, gdx_zeta);
	gradz(zeta, gdz_zeta);

	// calculate lap vort
	laplacian(zeta, lap_zeta);

	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			zeta_change[I(i,j)] = - u[I(i,j)] * gdx_zeta[I(i,j)] - w[I(i,j)] * gdz_zeta[I(i,j)] + nu * lap_zeta[I(i,j)] + total_torq[I(i,j)];
		}
	}
	
	// calculate theta part
	gradx(theta, gdx_theta);
	gradz(theta, gdz_theta);
	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			theta_change[I(i,j)] = - u[I(i,j)] * gdx_theta[I(i,j)] - w[I(i,j)] * gdz_theta[I(i,j)];
		}
	}
	
	// calculate rho part
	gradx(rho, gdx_rho);
	gradz(rho, gdz_rho);
	for(int j = 0; j < zpts; ++j) {
		for(int i = 0; i < xpts; ++i) {
			rho_change[I(i,j)] = - u[I(i,j)] * gdx_rho[I(i,j)] - w[I(i,j)] * gdz_rho[I(i,j)];
		}
	}
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
		}

		for(int rk = 0 ; rk < 4; ++rk) {

			switch (rk) {
				case 0:
					cal_change(zeta, theta, rho, zeta_k1, theta_k1, rho_k1);
					break;
				case 1:
					evolve(zeta,  zeta_k1,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k1, dt/2.0, theta_tmp);
					evolve(rho,   rho_k1,   dt/2.0, rho_tmp);
					
					cal_change(zeta_tmp, theta_tmp, rho_tmp, zeta_k2, theta_k2, rho_k2);

					break;
				case 2:
					evolve(zeta,  zeta_k2,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k2, dt/2.0, theta_tmp);
					evolve(rho,   rho_k2,   dt/2.0, rho_tmp);
					
					cal_change(zeta_tmp, theta_tmp, rho_tmp, zeta_k3, theta_k3, rho_k3);
					break;
				case 3:
					evolve(zeta,  zeta_k3,  dt/2.0, zeta_tmp);
					evolve(theta, theta_k3, dt/2.0, theta_tmp);
					evolve(rho,   rho_k3,   dt/2.0, rho_tmp);
					
					cal_change(zeta_tmp, theta_tmp, rho_tmp, zeta_k4, theta_k4, rho_k4);
					break;
			}
		}

		for(int j=0; j<zpts; ++j) {
			for(int i=0; i<xpts; ++i) {
				zeta[I(i,j)] += (zeta_k1[I(i, j)] + 2.0 * zeta_k2[I(i, j)] + 2.0 * zeta_k3[I(i, j)] + zeta_k4[I(i, j)]) / 6.0;
			}
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
	 * The input cannot be protected from damaging.
	 * Boundarz condition is specifz on the margin of input
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
		fprintf(stderr,"[%d] Err avg: %5f\n",iter, err);

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
/*
int main() {

	init();

	int err;

	auto evolution = [&] (ModelData& k, ModelData& evo_vort) {
		vort_tmp2.copzData(evo_vort);
		EnuData::each_index([&](xsize i, xsize j) {
			vort_tmp2(i,j) *= bc_filter(i,j);
		});

		// inverse vort to psi
		if((err = jacobi_relaxation(vort_tmp2, psi, tmp, dx, dz, max_err, max_iter)) != 0) {
			fprintf(stderr, "[Error code %d]: Relaxation cannot converge within [%d] iterations.\n", err, max_iter);
			return -1;
		}

		// calculate u,v
		gradz(psi, u, dz); u *= -1.0f;
		gradx(psi, v, dx);

		// calculate grad vort
		gradx(evo_vort, gradx_vort, dx);
		gradz(evo_vort, gradz_vort, dz);

		// calculate lap vort
		lap(evo_vort, lap_vort, dx, dz);
		
		EnuData::each_index([&](xsize i, xsize j){
			k(i,j) = -u(i,j) * gradx_vort(i,j) - v(i,j) * gradz_vort(i,j) + nu * lap_vort(i,j);	
		});
	
	};

	auto calculator = [&] (ModelData& output, ModelData& input, ModelData& k, T dt) {
		EnuData:: each_index([&](xsize i, xsize j) {
			output(i,j) = input(i,j) + k(i,j) * dt;		
		});
	};

	for(xsize i = 0; i < step; ++i) {
		cout << "step:" << i << endl;
		evolution(k1, vort); 	 calculator(vort_tmp, vort, k1, 0.5 * dt);
		evolution(k2, vort_tmp); calculator(vort_tmp, vort, k2, 0.5 * dt);
		evolution(k3, vort_tmp); calculator(vort_tmp, vort, k3, dt);
		evolution(k4, vort_tmp);

		EnuData::each_index([&](xsize i, xsize j){
			vort(i,j) +=  dt / 6.0 * (k1(i,j) + 2*k2(i,j) + 2*k3(i,j) + k4(i,j));
		});


		if(i % 10 == 0) {
			stringstream filename;
			filename << "vort_" << i << ".bin";
			output_data(vort, filename.str().c_str());
		}
		// calculate new vort
	}
//	cout << "start" << endl;
//	cout << "Jacobi Result: " << jacobi_relaxation(vort, psi, tmp, dx, dz, max_err, max_iter) << endl;
	//lap(psi,vort,dx,dz);

//	psi = vort;

//	cout << &(psi.getData(0)) << " compare " << (psi.getDataPointer()) << endl;

//	cout << "Seg?" <<endl;
//	psi(0,0) += 200;


//	output_data(vort, "init_vort.bin");
//	output_data(psi, "init_psi.bin");

	return 0;
}

*/
