#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        void test_quadrature(int order, bool border)
        { 
        	Quadrature quad = Quadrature::get_quadrature(order, border);
        	std::cout << quad.nb_points() << std::endl;
        	double sum = 0;
        	for ( int i = 0; i < quad.nb_points(); ++i) {
        		std::cout << quad.point(i).x << " " << quad.point(i).y << std::endl;
        		std::cout << quad.weight(i) << std::endl;
        		sum += quad.weight(i);
        	}
        	std::cout << sum << std::endl;
        }
        
        void test_mapping()
        {
        	Mesh mesh;
        	mesh.load("data/square.mesh"); //récupère le cas square
        	ElementMapping mapping(mesh, false, 4);
        	vertex point_ref; 
        	point_ref.x = 0.2; 
        	point_ref.y = 0.4;
        	
        	// Test de la méthode transform
        	std::cout << mapping.transform(point_ref).x << " " << mapping.transform(point_ref).y << '\n';
        	
        	// Test de la méthode jacobian_matrix
        	mapping.jacobian_matrix(point_ref).print();
        	
        	// Test de la méthode jacobian
        	std::cout << mapping.jacobian(point_ref) << '\n';
        }
        
        void test_ShapeFunc()
        {
        	vertex point_ref; 
        	point_ref.x = 0.1; 
        	point_ref.y = 0.6;
        	
        	// Cas d'un segment
        	ShapeFunctions shape_f_segment(1,1);
        	std::cout << shape_f_segment.nb_functions() << '\n';
        	
        	// Cas d'un triangle
        	ShapeFunctions shape_f_triangle(2,1);
        	std::cout << shape_f_triangle.nb_functions() << '\n';
        }
        
        void test_Ke()
        {
        }
        
        void test_src_term()
        {
        }

    }
}
