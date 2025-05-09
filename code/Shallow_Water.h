#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H
///////////////////////////////////////////////////////////////////////
//                  Shallow_Water.h
//                  MCG4139 Winter 2021
///////////////////////////////////////////////////////////////////////
#include<Eigen/Core>
#include<cmath>

///////////////////////////////////////////////////////////////////////
//                   Shallow_Water_Equations
///////////////////////////////////////////////////////////////////////
class Shallow_Water_Equations {
public:

  using Vector_type = Eigen::Vector3d;
  static constexpr int number_of_unknowns = 3;
  
  ////////////////////////////////////
  // get_primitive
  template<typename Vec_in>
  static Vector_type get_conserved(const Vec_in& U, const int& cell_id){
    if(cell_id == 200){
      return {U[0],
	    U[1]/U[0],
	    U[2]/U[0]};
    } else {
      Vector_type U_zero;
      U_zero.setZero();
      return U_zero;
    }
    
  }

  ////////////////////////////////////
  //  Flux x
  template<typename Vec_in>
  static Vector_type Fx(const Vec_in& U) {
    Vector_type F;
    F[0] = U[1];
    F[1] = U[1]*U[1]/U[0] + 0.5*g*U[0]*U[0];
    F[2] = U[1]*U[2]/U[0];
    return F;
  }

  ////////////////////////////////////
  //  Rotate
  template<typename Vec_in>
  static Vector_type rotate(const Vec_in& U, Node2D n_hat) {
    Vector_type U_rot;
    U_rot[0] =  U[0];
    U_rot[1] =  U[1]*n_hat.x() + U[2]*n_hat.y();
    U_rot[2] = -U[1]*n_hat.y() + U[2]*n_hat.x();
    return U_rot;
  }

  ////////////////////////////////////
  //  Rotate back
  template<typename Vec_in>
  static Vector_type rotate_back(const Vec_in& U, Node2D n_hat) {
    Vector_type U_rot;
    U_rot[0] = U[0];
    U_rot[1] = U[1]*n_hat.x() - U[2]*n_hat.y();
    U_rot[2] = U[1]*n_hat.y() + U[2]*n_hat.x();
    return U_rot;
  }

  ////////////////////////////////////
  //  reflect_x
  template<typename Vec_in>
  static Vector_type reflect_x(const Vec_in& U) {
    Vector_type U_ref = U;
    U_ref[1] *= -1.0;
    return U_ref;
  }

  ////////////////////////////////////
  //  inlet_x (also just reflection)
  template<typename Vec_in>
  static Vector_type inlet_x(const Vec_in& U) {
    Vector_type U_in;
    U_in[0] = U[0];
    U_in[1] = U[1];//0.0;//1.0e-3;
    U_in[2] = 0.0;   
    return U_in;
  }

  ////////////////////////////////////
  //  outlet_x (also just reflection)
  template<typename Vec_in>
  static Vector_type outlet_x(const Vec_in& U) {
    Vector_type U_out;
    U_out[0] = U[0];
    U_out[1] = U[1];
    U_out[2] = 0.0;
    return U_out;
  }

  ////////////////////////////////////
  //  no_slip_wall_x (just another reflection for shallow water eq.)
  template<typename Vec_in>
  static Vector_type no_slip_wall_x(const Vec_in& U) {
    Vector_type U_wall = U;
    const double ff = 1.0;
    U_wall[1]  *= -ff;
    return U_wall;
  }

  ////////////////////////////////////
  //       Source term
  template<typename Vec_in, typename EM_in>
  static Vector_type source(const Vec_in& U, const EM_in& em) {
    Vector_type S;
    const double source_multiplier = U[0]*magnetic_saturization*magnetic_fraction/rho;
    S[0] = 0.0;
    S[1] = source_multiplier*em[0];
    S[2] = source_multiplier*em[1];
    return S;
  }

  ////////////////////////////////////
  //  max_lambda_x
  template<typename Vec_in>
  static double max_lambda_x(const Vec_in& U) {
    return fabs(U[1]/U[0])+sqrt(U[0]*g);
  }

private:
  static constexpr double               tolerance = 1.0e-12;
  static constexpr double                       g = 9.81;
  static constexpr double                     rho = 1000;
  static constexpr double   magnetic_saturization = 375000.0;
  static constexpr double       magnetic_fraction = 0.1;
  static constexpr double permeability_free_space = (4*acos(-1.0))*0.0000001;
};


#endif //#ifndef SHALLOW_WATER_H
