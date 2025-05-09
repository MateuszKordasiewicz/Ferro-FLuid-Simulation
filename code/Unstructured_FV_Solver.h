#ifndef UNSTRUCTURED_FV_SOLVER_H
#define UNSTRUCTURED_FV_SOLVER_H
///////////////////////////////////////////////////////////////////////
//                  Unstructured_FV_Solver.h
//                  MCG4139 Winter 2021
///////////////////////////////////////////////////////////////////////
#include<cmath>
#include<string>
#include<functional>
#include<Eigen/Core>
#include"Unstructured_Mesh.h"
#include"EM_Background.h"

///////////////////////////////////////////////////////////////////////
//                    Local Lax-Friedrichs
///////////////////////////////////////////////////////////////////////
template<typename PDE_type, typename Vec_type>
auto Local_Lax_Friedrichs(const Vec_type& Ul, const Vec_type& Ur) {
  const auto Fl = PDE_type::Fx(Ul);
  const auto Fr = PDE_type::Fx(Ur);
  const auto max_lambda = std::max(PDE_type::max_lambda_x(Ul),PDE_type::max_lambda_x(Ur));
  typename PDE_type::Vector_type F = 0.5*(Fl+Fr-max_lambda*(Ur-Ul));
  return F;
}

///////////////////////////////////////////////////////////////////////
//             Unstructured Finite-Volume Scheme
///////////////////////////////////////////////////////////////////////
template<typename PDE_type>
class Unstructured_FV_Solver : public EM_Background {
public:

  
  using SolutionVector_type = typename PDE_type::Vector_type;
  using BackgroundVector_type = EM_Background::Vector_type;

  ///////////////////////////////////////////////////////////////////////
  //  Default Constructors, etc.
  Unstructured_FV_Solver() = default;
  Unstructured_FV_Solver(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver(Unstructured_FV_Solver&&) = default;
  Unstructured_FV_Solver& operator=(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver& operator=(Unstructured_FV_Solver&&) = default;

  ///////////////////////////////////////////////////////////////////////
  //  Constructor taking a filename and initial condition
  Unstructured_FV_Solver(const std::string& mesh_filename,
                         const std::function<SolutionVector_type(const Node2D&)>& ic);

  ///////////////////////////////////////////////////////////////////////
  //  Solution vector in cell "i"
  auto& U(int i) {
    return Global_U[i];
  }

  ///////////////////////////////////////////////////////////////////////
  //  Solution vector in cell "i"
  auto& S(int i) {
    return Global_S[i];
  }
  
  ///////////////////////////////////////////////////////////////////////
  //  dUdt in cell "i"
  auto& dUdt(int i) {
    return Global_dUdt[i];
  }

  ///////////////////////////////////////////////////////////////////////
  //  Background  vector in cell "i"
  auto& EM(int i) {
    return Global_EM[i];
  }

  auto& EM() {
    return Global_EM;
  }

  auto& centroids(int i){
    return Global_centroids[i];
  }

  ///////////////////////////////////////////////////////////////////////
  //  number of cells, nodes, and edges.
  auto number_of_cells() {return static_cast<int>(cells.size());}
  auto number_of_nodes() {return static_cast<int>(nodes.size());}
  auto number_of_edges() {return static_cast<int>(edges.size());}
  auto number_of_background_edges() {return static_cast<int>(background_edges.size());}

  ///////////////////////////////////////////////////////////////////////
  //  time march to time
  void time_march_to_time(double final_time, double CFL);

  ///////////////////////////////////////////////////////////////////////
  //  Make movie
  void make_movie(double final_time, double CFL, int num_frames, const std::string& file_location, const std::string& filename);

  ///////////////////////////////////////////////////////////////////////
  //  Output two vector max and fluid solution vector
  void Output_Vector_min_max(const std::string& file_location, const std::string& filename);
  void Output_Vector(const std::string& file_location, const std::string& filename);

  ///////////////////////////////////////////////////////////////////////
  //  write to VTK
  void write_to_vtk(const std::string& filename);

private:

  ///////////////////////////////////////////////////////////////////////
  //  Member variables
  double                                    h_vel;
  double                                    h_amp;
  double                                    e_amp;
  std::vector<Node2D>                       nodes;
  std::vector<Cell>                         cells;
  std::vector<Edge>                         edges;
  std::vector<Edge>              background_edges;
  std::vector<double>                       areas;
  std::vector<Node2D>            Global_centroids;
  std::vector<SolutionVector_type>       Global_U;
  std::vector<SolutionVector_type>       Global_S;
  std::vector<SolutionVector_type>    Global_dUdt;
  std::vector<BackgroundVector_type>    Global_EM;
  //EM_base_vector_type              Global_EM_base;
  double                                     time;

  ///////////////////////////////////////////////////////////////////////
  //  compute all areas
  void compute_areas() {
    areas.resize(number_of_cells());
    for(int i = 0; i < number_of_cells(); ++i) {
      const Cell& cell = cells[i];
      const Node2D& n0 = nodes[cell.node0];
      const Node2D& n1 = nodes[cell.node1];
      const Node2D& n2 = nodes[cell.node2];
      areas[i] = 0.5*fabs( n0.x()*n1.y()-n0.y()*n1.x()
                          +n1.x()*n2.y()-n1.y()*n2.x()
                          +n2.x()*n0.y()-n2.y()*n0.x());
    }
  }
};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
//  Constructor
template<typename PDE_type>
Unstructured_FV_Solver<PDE_type>::Unstructured_FV_Solver(const std::string& mesh_filename,
                                                         const std::function<SolutionVector_type(const Node2D&)>& ic) {

  ///////////////////////////////////////////////////////////////
  // Read gmsh file & compute mesh vectors for solution
  std::tie(nodes, cells) = read_gmsh_file(mesh_filename);
  std::tie(edges, background_edges) = compute_edges(nodes,cells);

  ///////////////////////////////////////////////////////////////
  // Diagnostics
  //print_edges(edges);
  //print_cells(cells);

  /////////////////////////////////////////////////////////////////////////////////
  // Resizing of vectors for our solver
  Global_U.resize(number_of_cells());
  Global_dUdt.resize(number_of_cells());
  Global_S.resize(number_of_cells());
  Global_EM.resize(number_of_cells());
  Global_centroids.resize(number_of_cells());
  ////////////////////////////////////////////////////////////////////////////////
  // Initialization
  time = 0.0;
  for(int i = 0; i < number_of_cells(); ++i) {
    centroids(i) = (nodes[cells[i].node0]+nodes[cells[i].node1]+nodes[cells[i].node2])/3.0;
    if(cells[i].id == 200){
      U(i) = ic(centroids(i)); 
    } else {
      SolutionVector_type U_zero;
      U_zero.setZero();
      U(i) = U_zero;
    }
  }
  EM() = this->EM_Background::solution(cells, Global_centroids);
  compute_areas();
}

///////////////////////////////////////////////////////////////////////
//  time march to time
template<typename PDE_type>
void Unstructured_FV_Solver<PDE_type>::time_march_to_time(double final_time, double CFL) {

  constexpr double tolerance = 1.0e-12;
  
  while(time < final_time - tolerance) {

    double dt = final_time-time;
    
    for(auto& entry : Global_dUdt) {
      entry.fill(0.0);
    }

    for(const auto& edge : edges) {

      double l;
      Node2D n_hat;
      std::tie(n_hat,l) = n_hat_and_length(nodes[edge.node0], nodes[edge.node1]);

      typename PDE_type::Vector_type Ul;
      typename PDE_type::Vector_type Ur;
      typename PDE_type::Vector_type Ul_rot;
      typename PDE_type::Vector_type Ur_rot;
      
      if(edge.cell_l != -100 && edge.cell_l != -101 && edge.cell_l != -102 && edge.cell_l != -105) {
        Ul = U(edge.cell_l);
        Ul_rot = PDE_type::rotate(Ul, n_hat);
      }

      if(edge.cell_r != -100 && edge.cell_r != -101 && edge.cell_r != -102 && edge.cell_r != -105) {
        Ur = U(edge.cell_r);
        Ur_rot = PDE_type::rotate(Ur, n_hat);
      }

      if(edge.cell_l == -100) {
        Ul_rot = PDE_type::inlet_x(Ur_rot);
	//Ul_rot = PDE_type::reflect_x(Ur_rot);
      }

      if(edge.cell_r == -100) {
        Ur_rot = PDE_type::inlet_x(Ul_rot);
	//Ur_rot = PDE_type::reflect_x(Ul_rot);
      }

      if(edge.cell_l == -101) {
        Ul_rot = PDE_type::outlet_x(Ur_rot);
	//Ul_rot = PDE_type::reflect_x(Ur_rot);
      }

      if(edge.cell_r == -101) {
        Ur_rot = PDE_type::outlet_x(Ul_rot);
	//Ur_rot = PDE_type::reflect_x(Ul_rot);
      }

      if(edge.cell_l == -102) {
        Ul_rot = PDE_type::no_slip_wall_x(Ur_rot);
      }

      if(edge.cell_r == -102) {
        Ur_rot = PDE_type::no_slip_wall_x(Ul_rot);
      }
      
      auto F_rot = Local_Lax_Friedrichs<PDE_type>(Ul_rot, Ur_rot);
      auto F = PDE_type::rotate_back(F_rot, n_hat);
      
 
      if(edge.cell_l != -100 && edge.cell_l != -101 && edge.cell_l != -102 && edge.cell_l != -105) {
	dUdt(edge.cell_l) -= F*l/areas[edge.cell_l];
	double max_lambda = PDE_type::max_lambda_x(Ul);
	dt = std::min(dt, CFL*sqrt(areas[edge.cell_l])/max_lambda);
      }

      if(edge.cell_r != -100 && edge.cell_r != -101 && edge.cell_r != -102 && edge.cell_r != -105) {
	dUdt(edge.cell_r) += F*l/areas[edge.cell_r];
	double max_lambda = PDE_type::max_lambda_x(Ur);
	dt = std::min(dt, CFL*sqrt(areas[edge.cell_r])/max_lambda);
      }

    } //end loop over edges
    
    for(int i = 0; i < number_of_cells(); ++i) {
      if(cells[i].id == 200){
        typename PDE_type::Vector_type source = PDE_type::source(U(i), EM(i));
	U(i) += dt*(dUdt(i) + source);
      }
    }
    
    time += dt;
    std::cout << "Time = " << time << "    dt = " << dt << '\n';
  }
}


///////////////////////////////////////////////////////////////////////
//  Make movie
template<typename PDE_type>
void Unstructured_FV_Solver<PDE_type>::make_movie(double final_time,
                                                  double CFL,
                                                  int num_frames,
                                                  const std::string& file_location,
						  const std::string& filename) {

  const auto dt = final_time/static_cast<double>(num_frames-1);
  
  const std::string filename_base = file_location + "movie/movie" + filename + "_";

  auto get_name = [&filename_base] (int i) {
    std::stringstream ss;
    ss << filename_base << std::setfill('0') << std::setw(10) << i << ".vtk";
    return ss.str();
  };

  write_to_vtk(get_name(0));

  for(int i = 1; i < num_frames; ++i) {
    std::cout << "Time Marching to time = " << dt*static_cast<double>(i) << '\n';
    time_march_to_time(dt*static_cast<double>(i), CFL);
    write_to_vtk(get_name(i));
  }
  Output_Vector_min_max(file_location, filename);
  Output_Vector(file_location, filename);
}

///////////////////////////////////////////////////////////////////////
//           Output cell_types to .dat file (diagnostic)
///////////////////////////////////////////////////////////////////////
auto print_cell_types(const std::vector<int>& cell_types) {
  
  std::string filename = "cell_types.dat";
  std::cout << "Writing output to " << filename << ".\n";
  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  
  fout.precision(16);
  fout << std::setw(30) << "cell #"
       << std::setw(30) << "type" << "\n";
  
  for(int i = 0; i < static_cast<int>(cell_types.size()); ++i){
    fout << std::setw(30) << i
	 << std::setw(30) << cell_types[i] << "\n";  
  }
}

///////////////////////////////////////////////////////////////////////
//           Output Solution vector min and max to .dat
///////////////////////////////////////////////////////////////////////
template<typename PDE_type>
void Unstructured_FV_Solver<PDE_type>::Output_Vector_min_max(const std::string& file_location, const std::string& filename) {
  
  const std::string filename_base = file_location + "min_max/vector_min_max" + filename + ".dat";

  typename PDE_type::Vector_type maxValues = PDE_type::Vector_type::Zero();
  typename PDE_type::Vector_type minValues = PDE_type::Vector_type::Zero();

  for(int i = 0; i < number_of_cells(); ++i){
    if(cells[i].id == 200){
      const auto W = PDE_type::get_conserved(U(i), cells[i].id);
    
      maxValues[0] = std::max(maxValues[0], W[0]);
      maxValues[1] = std::max(maxValues[1], W[1]);
      maxValues[2] = std::max(maxValues[2], W[2]);
    
      minValues[0] = std::min(maxValues[0], W[0]);
      minValues[1] = std::min(maxValues[1], W[1]);
      minValues[2] = std::min(maxValues[2], W[2]);
    }
  }
  
  std::cout << "Writing min and max to " << filename_base << "\n";
  std::ofstream fout(filename_base);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename_base);
  }
  
  fout.precision(16);
  fout << std::setw(30) << "Height_min"
       << std::setw(30) << "Height_max"
       << std::setw(30) << "X_velocity_min"
       << std::setw(30) << "X_velocity_max"
       << std::setw(30) << "Y_velocity_min"
       << std::setw(30) << "Y_velocity_max" << "\n";
  
  fout << std::setw(30) << minValues[0]
       << std::setw(30) << maxValues[0]
       << std::setw(30) << minValues[1]
       << std::setw(30) << maxValues[1]
       << std::setw(30) << minValues[2]
       << std::setw(30) << maxValues[2] << "\n";  
  
}
    

///////////////////////////////////////////////////////////////////////
//           Output Solution vector  to .dat
///////////////////////////////////////////////////////////////////////
template<typename PDE_type>
void Unstructured_FV_Solver<PDE_type>::Output_Vector(const std::string& file_location, const std::string& filename) {
  
  const std::string filename_base = file_location + "vector/vector" + filename + ".dat";

  std::cout << "Writing vector to " << filename_base << "\n";
  std::ofstream fout(filename_base);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename_base);
  }
  
  fout.precision(16);
  
  fout << std::setw(30) << "X"
       << std::setw(30) << "Y"
       << std::setw(30) << "Height"
       << std::setw(30) << "X_velocity"
       << std::setw(30) << "Y_velocity" << "\n";;
  
  for(int i = 0; i < number_of_cells(); ++i){
    if(cells[i].id == 200){
      const auto W = PDE_type::get_conserved(U(i), cells[i].id);
      
      fout << std::setw(30) << centroids(i).x()
	   << std::setw(30) << centroids(i).y()
	   << std::setw(30) << W[0]
	   << std::setw(30) << W[1]
	   << std::setw(30) << W[2] << "\n";  
    }
  }
}

///////////////////////////////////////////////////////////////////////
//  write to VTK
template<typename PDE_type>
void Unstructured_FV_Solver<PDE_type>::write_to_vtk(const std::string& filename) {

  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("error opening vtk file.");
  }

  fout << "# vtk DataFile Version 2.0\n"
       << "Unstructured Solver\n"
       << "ASCII\n"
       << "DATASET UNSTRUCTURED_GRID\n"
       << "POINTS " << number_of_nodes() << " double\n";

  for(const auto& node : nodes) {
    fout << node.x() << " " << node.y() << " 0\n";
  }

  fout << "\nCELLS " << number_of_cells() << " " << 4*number_of_cells() << '\n';
  for(const auto& cell : cells) {
    fout << "3 " << cell.node0 << " " << cell.node1 << " " << cell.node2 << '\n';
  }

  fout << "\nCELL_TYPES " << number_of_cells() << '\n';
  for(int i = 0; i < number_of_cells(); ++i) {
    fout << "5\n";
  }

  fout << "\nCELL_DATA " << number_of_cells()  << '\n'
       << "SCALARS Height double 1\n"
       << "LOOKUP_TABLE default\n";

  for(int i = 0; i < number_of_cells(); ++i) {
    const auto W = PDE_type::get_conserved(U(i), cells[i].id);
    fout << W[0] << '\n';
  }

  fout << "VECTORS Velocity double\n";
  
  for(int i = 0; i < number_of_cells(); ++i) {
    const auto W = PDE_type::get_conserved(U(i), cells[i].id);
    fout << W[1] << " " << W[2] << " 0.0" << '\n';
  }
  /////////////////////////////////////////////////////////////
  // Extra outputs (1.5x - 2x the file size)
  //fout << "VECTORS Magnetic_field double\n";
  
  //for(int i = 0; i < number_of_cells(); ++i) {
  //  fout << EM(i)[0] << " " << EM(i)[1] << " 0.0" << '\n';
  //}

  //fout << "SCALARS Cell_type double 1\n"
  //     << "LOOKUP_TABLE default\n";

  //std::vector<int> cell_types(number_of_cells(), 0);
  //
  //for(int i = 0; i < number_of_edges(); ++i) {
  //  const Edge edge = edges[i];
  //  if(edge.cell_l == -100){
  //    cell_types[edge.cell_r] = -100;
  // } 
  //  if(edge.cell_l == -101){
  //    cell_types[edge.cell_r] = -200;
  //  }
  //  if(edge.cell_l == -102){
  //    cell_types[edge.cell_r] = -300;
  //  }
  //  if(edge.cell_r == -100){
  //    cell_types[edge.cell_l] = 100;
  //  }
  //  if(edge.cell_r == -101){
  //    cell_types[edge.cell_l] = 200;
  //  }
  //  if(edge.cell_r == -102){
  //    cell_types[edge.cell_l] = 300;
  // }
  //}

  //for(int i = 0; i < number_of_background_edges(); ++i) {
  //  const Edge edge = background_edges[i];
  //  if(edge.cell_l == -105){
  //    cell_types[edge.cell_r] = -500;
  //  }
  //  if(edge.cell_r == -105){
  //    cell_types[edge.cell_l] = 500;
  //  }
  //}

  //for(int i = 0; i < number_of_cells(); ++i) {
  //  fout << cell_types[i] << '\n';
  //}

  //fout << "SCALARS Cell_domain double 1\n"
  //     << "LOOKUP_TABLE default\n";

  //for(int i = 0; i < number_of_cells(); ++i) {
  //  fout << cells[i].id << '\n';
  //}
}







#endif //#ifndef UNSTRUCTURED_FV_SOLVER_H
