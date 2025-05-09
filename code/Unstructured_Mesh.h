#ifndef UNSTRUCTURED_MESH_H
#define UNSTRUCTURED_MESH_H
///////////////////////////////////////////////////////////////////////
//                  Unstructured_Mesh.h
//                  MCG4139 Winter 2021
///////////////////////////////////////////////////////////////////////
#include<vector>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<stdexcept>
#include<sstream>
#include<Eigen/Core>
///////////////////////////////////////////////////////////////////////
//                Some useful types (useful everywhere)
///////////////////////////////////////////////////////////////////////
struct Cell {
  int node0;
  int node1;
  int node2;
  int line01_type;
  int line12_type;
  int line20_type;
  int id;
};
///////////////////////////////////////////////////////////////////////
struct Edge {
  int node0;
  int node1;
  int cell_l;
  int cell_r;
};
///////////////////////////////////////////////////////////////////////
struct Boundary {
  int node0;
  int node1;
  int type;
};
///////////////////////////////////////////////////////////////////////
using Node2D = Eigen::Vector2d;

///////////////////////////////////////////////////////////////////////
//   n_hat_and_length (usefull function [anywhere given two nodes])
///////////////////////////////////////////////////////////////////////
auto n_hat_and_length(const Node2D& n0, const Node2D& n1) {
  
  auto length = (n1-n0).norm();
  Node2D n_hat;
  n_hat.x() =  (n1.y()-n0.y())/length;
  n_hat.y() = -(n1.x()-n0.x())/length;
  return std::make_pair(n_hat, length);
}

///////////////////////////////////////////////////////////////////////
//           Output BC to .dat file (diagnostic)
///////////////////////////////////////////////////////////////////////
auto print_BC(const std::vector<Boundary>& BC) {
  
  std::string filename = "BC.dat";
  std::cout << "Writing output to " << filename << ".\n";
  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  
  fout.precision(16);
  fout << std::setw(30) << "node0"
       << std::setw(30) << "node1"
       << std::setw(30) << "type" << "\n";
  
  for(int i = 0; i < static_cast<int>(BC.size()); ++i){
    Boundary b = BC[i];
    fout << std::setw(30) << b.node0
	 << std::setw(30) << b.node1
	 << std::setw(30) << b.type << "\n";  
  }
}

///////////////////////////////////////////////////////////////////////
//        Output cells to .dat file (diagnostic)
//////////////////////////////////////////////////////////////////////
auto print_cells(const std::vector<Cell>& cells) {

  std::string filename = "cells.dat";
  std::cout << "Writing output to " << filename << ".\n";
  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  fout.precision(16);
  
  fout << std::setw(30) << "node0"
       << std::setw(30) << "node1"
       << std::setw(30) << "node2"
       << std::setw(30) << "line01_type"
       << std::setw(30) << "line12_type"
       << std::setw(30) << "line20_type" << "\n";  
  
  for(int i = 0; i < static_cast<int>(cells.size()); ++i){
    Cell cell = cells[i];
    fout << std::setw(30) << cell.node0
	 << std::setw(30) << cell.node1
	 << std::setw(30) << cell.node2
	 << std::setw(30) << cell.line01_type
	 << std::setw(30) << cell.line12_type
	 << std::setw(30) << cell.line20_type << "\n";  
  }
}

///////////////////////////////////////////////////////////////////////
//                 Read gmsh file (for solver)
///////////////////////////////////////////////////////////////////////
auto read_gmsh_file(const std::string& filename) {
  
  std::ifstream fin(filename);
  if(!fin) {
    throw std::runtime_error("Cannot open file: " + filename + ".");
  }
  ///////////////////////////////////////////////////////
  // Lambda function to consume lines that are expected
  // but should be ignored
  auto expect_line = [filename] (std::ifstream& fin, const std::string& expected) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    if(s != expected) {
      throw std::runtime_error("Error reading file: " + filename + ".\n" +
                               "Expected \"" + expected + "\", but got \"" + s +
                               "\".");
    }
  };
  std::cout << "\n\nReading gmsh file: " << filename << '\n';
  
  ///////////////////////////////////////////////////////
  // Read file
  expect_line(fin, "$MeshFormat");
  expect_line(fin, "2.2 0 8");
  expect_line(fin, "$EndMeshFormat");
  expect_line(fin, "$Nodes");

  int number_of_nodes;
  fin >> number_of_nodes;

  std::vector<Node2D> nodes(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i) {
    int    dummy_index;
    double dummy_z_coordinate;
    fin >> dummy_index
        >> nodes[i].x()
        >> nodes[i].y()
        >> dummy_z_coordinate;
    if(dummy_index != i+1) {
      throw std::runtime_error("Error with node index.");
    }
    if(fabs(dummy_z_coordinate) > 1.0e-12) {
      throw std::runtime_error("Error, node has z component.");
    }
  }

  expect_line(fin, "$EndNodes");
  expect_line(fin, "$Elements");

  int number_of_elements; //not all will be cells
  fin >> number_of_elements;

  std::vector<Cell> cells;
  std::vector<Boundary> BC;
  for(int i = 0; i < number_of_elements; ++i) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    std::istringstream ss(s);

    int element_num;
    ss >> element_num;
    if(element_num - 1 != i) {
      throw std::runtime_error("Error reading element number.");
    }

    int element_type;
    ss >> element_type;
    
    if(element_type == 1){
      Boundary tmp_b;
      int dummy;
      ss >> dummy >> tmp_b.type >> dummy >> tmp_b.node0 >> tmp_b.node1;
      tmp_b.node0 -= 1;
      tmp_b.node1 -= 1;
      BC.push_back(tmp_b);
    }

    if(element_type == 2) {
      //triangular cell
      Cell c;
      int dummy;
      int cell_id;
      ss >> dummy >> cell_id >> dummy >> c.node0 >> c.node1 >> c.node2; 
      c.node0 -= 1;
      c.node1 -= 1;
      c.node2 -= 1;
      c.id = cell_id;

      auto assign_bc = [&c](int node_a, int node_b, int &bc_value_a, const std::vector<Boundary> &bc_data) {
	auto it = std::find_if(bc_data.begin(), bc_data.end(), [node_a, node_b](const Boundary &bc){
	  if(bc.node0 == node_a && bc.node1 == node_b) return true;
	  if(bc.node0 == node_b && bc.node1 == node_a) return true;
	  return false;
	});
	if(it != bc_data.end()){
	  bc_value_a = -it->type;
	}else {
	  bc_value_a = -1;
	}
      };
      
      assign_bc(c.node0, c.node1, c.line01_type, BC);
      assign_bc(c.node1, c.node2, c.line12_type, BC);
      assign_bc(c.node2, c.node0, c.line20_type, BC);
	   
      cells.push_back(c);
    }
  }
  ////////////////////////////
  // Diagnostic
  //print_BC(BC);
  //print_cells(cells);
  
  expect_line(fin, "$EndElements");

  std::cout << "done.\n";

  return std::make_pair(nodes, cells);

}

///////////////////////////////////////////////////////////////////////
//         Compute Edges (for read_gmsh -> for solver)
///////////////////////////////////////////////////////////////////////
auto compute_edges(const std::vector<Node2D>& nodes, const std::vector<Cell>& cells) {

  std::cout << "Computing Edges.\n";
  double target = 0.01;

  std::vector<Edge> edges;
  std::vector<Edge> background_edges;
  
  for(int i = 0; i < static_cast<int>(cells.size()); ++i) {

    const auto& cell = cells[i];
    
    /////////////////////////////////////////////////////////////////////////
    //   fluid domain edges
    if(cell.id == 200){
      auto e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
	if(cell.node0 == edge.node0 && cell.node1 == edge.node1) return true;
	if(cell.node0 == edge.node1 && cell.node1 == edge.node0) return true;
	return false;
      });

      if(e != edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	edges.push_back({cell.node0, cell.node1, i, cell.line01_type});
      }

      e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
	if(cell.node1 == edge.node0 && cell.node2 == edge.node1) return true;
	if(cell.node1 == edge.node1 && cell.node2 == edge.node0) return true;
	return false;
      });

      if(e != edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	edges.push_back({cell.node1, cell.node2, i, cell.line12_type});
      }

      e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
	if(cell.node2 == edge.node0 && cell.node0 == edge.node1) return true;
	if(cell.node2 == edge.node1 && cell.node0 == edge.node0) return true;
	return false;
      });

      if(e != edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	edges.push_back({cell.node2, cell.node0, i, cell.line20_type});
      }
    } else if(cell.id == 100){
      /////////////////////////////////////////////////////////////////////////
      //    background edges
      auto e = std::find_if(background_edges.begin(), background_edges.end(), [&cell] (const Edge& edge) {
	if(cell.node0 == edge.node0 && cell.node1 == edge.node1) return true;
	if(cell.node0 == edge.node1 && cell.node1 == edge.node0) return true;
	return false;
      });

      if(e != background_edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	background_edges.push_back({cell.node0, cell.node1, i, cell.line01_type});
      }

      e = std::find_if(background_edges.begin(), background_edges.end(), [&cell] (const Edge& edge) {
	if(cell.node1 == edge.node0 && cell.node2 == edge.node1) return true;
	if(cell.node1 == edge.node1 && cell.node2 == edge.node0) return true;
	return false;
      });

      if(e != background_edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	background_edges.push_back({cell.node1, cell.node2, i, cell.line12_type});
      }

      e = std::find_if(background_edges.begin(), background_edges.end(), [&cell] (const Edge& edge) {
	if(cell.node2 == edge.node0 && cell.node0 == edge.node1) return true;
	if(cell.node2 == edge.node1 && cell.node0 == edge.node0) return true;
	return false;
      });

      if(e != background_edges.end()) { //it was found
	e->cell_r = i;
      } else { // new edge
	background_edges.push_back({cell.node2, cell.node0, i, cell.line20_type});
      }

      if(static_cast<double>(i+1)/static_cast<double>(cells.size()) >= target) {
	std::cout << i+1 << "/" << cells.size() << " cells processed.\n";
	std::cout.flush();
	target += 0.01;
      }
    }
  }

  std::cout << "All cells processed.\n";

  //make sure all edges are properly aligned
  // (that the unit normal points from cell_l to cell_r)
  for(auto& edge: edges) {
      using std::swap;

      const auto edge_centroid = 0.5*(nodes[edge.node0]+nodes[edge.node1]);

      const auto& cell_left = cells[edge.cell_l];  //cell_left is never -1
      const auto left_cell_centroid = ( nodes[cell_left.node0]
                                       +nodes[cell_left.node1]
                                       +nodes[cell_left.node2])/3.0;

      const auto vec1 = left_cell_centroid - edge_centroid;

      auto n_hat_l = n_hat_and_length(nodes[edge.node0],nodes[edge.node1]);
      const auto& n_hat = n_hat_l.first;

      if(n_hat.dot(vec1) > 0.0) swap(edge.cell_l, edge.cell_r);
    
  }


  return std::make_pair(edges, background_edges);

}

///////////////////////////////////////////////////////////////////////
//        Output edges to .dat file (diagnostic)
//////////////////////////////////////////////////////////////////////
auto print_edges(const std::vector<Edge>& edges) {

  std::string filename = "edges.dat";
  std::cout << "Writing output to " << filename << ".\n";
  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  fout.precision(16);

  
  fout << std::setw(30) << "node0"
       << std::setw(30) << "node1"
       << std::setw(30) << "cell_l"
       << std::setw(30) << "cell_r" << "\n";  
  
  for(int i = 0; i < static_cast<int>(edges.size()); ++i){
    Edge edge = edges[i];
    fout << std::setw(30) << edge.node0
	 << std::setw(30) << edge.node1
	 << std::setw(30) << edge.cell_l
	 << std::setw(30) << edge.cell_r << "\n";
  }
}


#endif //o#ifndef UNSTRUCTURED_MESH_H
