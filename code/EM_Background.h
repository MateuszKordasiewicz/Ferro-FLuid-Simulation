#ifndef EM_BACKGROUND_H
#define EM_BACKGROUND_H
///////////////////////////////////////////////////////////////////////
//                  Shallow_Water.h
//                  MCG4139 Winter 2021
///////////////////////////////////////////////////////////////////////
#include<vector>
#include<Eigen/Core>
#include"Unstructured_Mesh.h"

struct wire{
  std::vector<int> cells;
  Node2D centroid;

  void set_zero(){
    cells.clear();
    centroid.setZero();
  }
};
  
///////////////////////////////////////////////////////////////////////
//                   Electric Magnetic Background
///////////////////////////////////////////////////////////////////////
class EM_Background {
public:
     
  ///////////////////////////////////////////////////////////////////////
  //  Default Constructors, etc.
  EM_Background() = default;
  EM_Background(const EM_Background&) = default;
  EM_Background(EM_Background&&) = default;
  EM_Background& operator=(const EM_Background&) = default;
  EM_Background& operator=(EM_Background&&) = default;
	
  using Vector_type = Eigen::Vector3d;
  static constexpr int number_of_unknowns = 3;

  ////////////////////////////////////
  //  Solution Initialization

  std::vector<Vector_type> solution(const std::vector<Cell>& edges, const std::vector<Node2D>& centroids);
  
private:  
  static constexpr double                 current = 1.0e2;       // Magnitude of the current
  static constexpr int         magnetic_direction = 1;       // Which way do we want the field lines to go - follows right hand rule - positive wires in +y will have current (*) if this is positive
  const double                      wire_diameter = 0.01;       // Needs to be changed based on wire diameter from the mesh
  static constexpr double permeability_free_space = (4*acos(-1.0))*0.0000001;
};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//          Generate magnetic field within domain based on the wires.
////////////////////////////////////////////////////////////////////////
std::vector<EM_Background::Vector_type> EM_Background::solution(const std::vector<Cell>& cells, const std::vector<Node2D>& centroids) {
  
  std::vector<Vector_type> magnetic_field;
  std::vector<wire> wires;
  std::vector<int> wire_cells;
  
  ////////////////////////////////////////////////////////////////////////
  //              Find cells bordering a wire
  for(int i = 0; i < static_cast<int>(cells.size()); ++i) {
    if(cells[i].id == 300) {
      wire_cells.push_back(i);
    }
  }
  
  ////////////////////////////////////////////////////////////////////////
  //            Function declerations
  auto distance = [](const Node2D& node0, const Node2D& node){
    const double x = (node0.x()-node.x());
    const double y = (node0.y()-node.y());
    const double distance = sqrt(x*x+y*y);
    //std::cout << "distance = " << distance << "\n" << std::flush;
    return distance;
  };

  auto find_centroid = [&centroids](const std::vector<int>& cells){
    Node2D wire_centroid;
    const int precision = 1; 
    wire_centroid.setZero();
    for(const auto& cell : cells){
      wire_centroid.x() += centroids[cell].x();
      wire_centroid.y() += centroids[cell].y();
    }
    wire_centroid.x() = std::round(wire_centroid.x()*std::pow(10, precision))/(cells.size()*std::pow(10, precision));
    wire_centroid.y() = std::round(wire_centroid.y()*std::pow(10, precision))/(cells.size()*std::pow(10, precision));
    return wire_centroid;
  };

  auto set_magnetic_field = [&wires, &magnetic_field, &centroids](){
    const auto PI = acos(-1.0);
    magnetic_field.resize(static_cast<int>(centroids.size()));
    for(int i = 0; i < static_cast<int>(centroids.size()); ++i) {
      for(auto& wire : wires){
	Vector_type local_field;
	const double radius = (wire.centroid - centroids[i]).norm();
	const auto cell_direction = (wire.centroid - centroids[i])/radius;
	auto field_direction = Eigen::Vector2d(cell_direction.y(), -cell_direction.x());
	if(magnetic_direction < 0){
	  field_direction.x() *= -1.0;
	  field_direction.y() *= -1.0;
	}
	if(centroids[i].y() < 0.0){
	  field_direction.x() *= -1.0;
	  field_direction.y() *= -1.0;
	}
	const double magnitude = permeability_free_space*current/(2*PI*radius);
	local_field[0] = magnitude*field_direction.x();
	local_field[1] = magnitude*field_direction.y();
	local_field[2] = 0.0;
	magnetic_field[i] += local_field; 
      }
    }
  };
  
  ////////////////////////////////////////////////////////////////////////
  //          Seperating wires and finding centroids
  wire new_wire;
  int node0;
  
  for(int i = 0; i < static_cast<int>(wire_cells.size()); ++i){
    
    if(i == 0){
      
      new_wire.cells.push_back(wire_cells[i]);
      wires.push_back(new_wire);
    }  
 
    if(i != 0){
      
      bool added_new_wire = false;
      
      for(auto& wire : wires){
      
	node0 = wire.cells[0];
       
	if(distance(centroids[node0], centroids[wire_cells[i]]) <= wire_diameter){
	  
	  wire.cells.push_back(wire_cells[i]);
	  added_new_wire = true;
	  break;
	}
      }
      
      if(!added_new_wire){
	
	new_wire.set_zero();
	new_wire.cells.push_back(wire_cells[i]);
	wires.push_back(new_wire);
      }
    }
  }
  
  for(auto& wire : wires){
    wire.centroid = find_centroid(wire.cells);
  }
  ////////////////////////////////////////////////////////////////////////
  //         Print the centroids (diagnostic)
  std::cout << "\n" << "# of wires: " << static_cast<int>(wires.size()) << "\n" << std::flush;
  //std::cout << "\nWire centroids:\n" << std::flush;
  //for(const auto& wire : wires){ 
  //  std::cout << "wire centroid x = " << wire.centroid.x() << " , y = " << wire.centroid.y() << "\n" << std::flush;
  //}
  std::cout << "\n" << std::flush;
  
  ////////////////////////////////////////////////////////////////////////
  //          Setup and return of the magnetic field
  set_magnetic_field();
  return magnetic_field;
}

#endif //#ifndef EM_BACKGROUND_H
