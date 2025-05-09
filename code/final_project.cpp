#include"Unstructured_FV_Solver.h"
#include"Unstructured_Mesh.h"
#include"Shallow_Water.h"


struct Mesh_Study {
  std::vector<std::string>             meshes;
  std::string                 location_output;
  std::vector<std::string>  output_file_names;
};

struct Time_Study {
  std::string                           mesh;
  std::vector<double>                   CFLs;
  std::string                location_output;
  std::vector<std::string> output_file_names;
};

///////////////////////////////////////////////////////////////////////
//                         Main
///////////////////////////////////////////////////////////////////////
int main() {
  
  ///////////////////////////////////////////////////////
  // Output console to .dat file for diagnostics
  //std::freopen("diagnostic.dat", "w", stdout);
  
  ///////////////////////////////////////////////////////
  // Fun time
  std::cout << " --------------------------------------------------------\n"
               " |                                                      |\n"
               " |          unstructured shallow-water solver           |\n"
               " |           ferro-fluid in a magnetic field            |\n"
               " |                                                      |\n"
               " --------------------------------------------------------\n\n";

  auto ics = [] (const Node2D& n) {

    Eigen::Vector3d U;
    U[1] = 0.0;
    U[2] = 0.0;
    if(n.squaredNorm() < 0.05) {
      U[0] = 1.5;
    } else {
      U[0] = 1.0;
    }
    return U;
  };

  auto icc = [] (const Node2D& n) {

    Eigen::Vector3d U;
    U[0] = 1.0;
    U[1] = 0.0;
    U[2] = 0.0;
    return U;
  };

  const double      final_time =         0.5;
  double                   CFL =         0.25;
  const int           n_frames =          401;
  std::string             mesh = "mesh_study/meshes/domain_0.0005.msh";//"mesh_study/meshes/domain_0.005.msh";
  std::string  location_output = "output/";//"mesh_study/";
  std::string output_file_name = "_CL0.0005_0.5s";//"_CharLen_0.005";

  const bool mesh_study_on = false;  // Careful make sure you have space (or reduce number of frames or mesh complexity)
  Mesh_Study mesh_study;
  mesh_study.meshes = {"mesh_study/meshes/domain_0.05.msh", "mesh_study/meshes/domain_0.025.msh", "mesh_study/meshes/domain_0.0175.msh", "mesh_study/meshes/domain_0.01.msh",
  		       "mesh_study/meshes/domain_0.0075.msh", "mesh_study/meshes/domain_0.005.msh", "mesh_study/meshes/domain_0.0025.msh", "mesh_study/meshes/domain_0.00175.msh",
  		       "mesh_study/meshes/domain_0.001.msh", "mesh_study/meshes/domain_0.00075.msh"};            // ~8 GiB
  //mesh_study.meshes = {"mesh_study/meshes/domain_0.0005.msh", "mesh_study/meshes/domain_0.00025.msh"}; // ~47 GiB
  //mesh_study.meshes = {"mesh_study/meshes/domain_0.000175.msh"};                                       // ~82 GiB
  
  
  mesh_study.location_output = "mesh_study/";
  mesh_study.output_file_names = {"_CharLen_0.05", "_CharLen_0.025","_CharLen_0.0175", "_CharLen_0.01",
  				  "_CharLen_0.0075", "_CharLen_0.005", "_CharLen_0.0025", "_CharLen_0.00175",
  				  "_CharLen_0.001", "_CharLen_0.00075"};
  //mesh_study.output_file_names = {"_CharLen_0.0005", "_CharLen_0.00025"};
  //mesh_study.output_file_names = {"_CharLen_0.000175"};

  
  const bool time_study_on = false; // Careful make sure you have space (or reduce number of frames or mesh complexity)
  Time_Study time_study;
  time_study.mesh = "time_study/domain.msh";
  time_study.CFLs = {1.01, 0.99, 0.9, 0.75, 0.5, 0.25};
  time_study.location_output = "time_study/";
  time_study.output_file_names = {"_CFL_1.01", "_CFL_0.99", "_CFL_0.9", "_CFL_0.75", "_CFL_0.5", "_CFL_0.25"}; // ~89 GiB 

  
  if(mesh_study_on){
    
    for(int i = 0; i < static_cast<int>(mesh_study.meshes.size()); ++i){
      
      mesh = mesh_study.meshes[i];
      
      auto solver = Unstructured_FV_Solver<Shallow_Water_Equations>(mesh, icc);
      solver.make_movie(final_time, CFL, n_frames, mesh_study.location_output, mesh_study.output_file_names[i]);
    }
  } else if(time_study_on){
    
    for(int i = 0; i < static_cast<int>(time_study.CFLs.size()); ++i){
      
      CFL = time_study.CFLs[i];
      mesh = time_study.mesh;
      
      auto solver = Unstructured_FV_Solver<Shallow_Water_Equations>(mesh, icc);
      solver.make_movie(final_time, CFL, n_frames, time_study.location_output, time_study.output_file_names[i]);
    }
  } else {
    
    auto solver = Unstructured_FV_Solver<Shallow_Water_Equations>(mesh, icc);
    solver.make_movie(final_time, CFL, n_frames, location_output, output_file_name);
  }
  
  return 0;
}

