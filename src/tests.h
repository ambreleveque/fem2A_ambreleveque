#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"
#include "simu.h"

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
        	std::cout << "Somme des poids : " << sum << std::endl;
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
        	point_ref.y = 0.4;
        	
        	// Cas d'un segment
        	ShapeFunctions shape_f_segment(1,1);
        	std::cout << shape_f_segment.nb_functions() << '\n';
        	std::cout << "Shape function d'indice 0 : " << shape_f_segment.evaluate(0,point_ref) << '\n';
        	std::cout << "Shape function d'indice 1 : " << shape_f_segment.evaluate(1,point_ref) << '\n';
        	std::cout << "Vecteur gradient de la shape function d'indice 0 : " << shape_f_segment.evaluate_grad(0,point_ref).x << " " << shape_f_segment.evaluate_grad(0,point_ref).y << '\n';
        	std::cout << "Vecteur gradient de la shape function d'indice 1 : " << shape_f_segment.evaluate_grad(1,point_ref).x << " " << shape_f_segment.evaluate_grad(1,point_ref).y << '\n';
        	
        	// Cas d'un triangle
        	ShapeFunctions shape_f_triangle(2,1);
        	std::cout << shape_f_triangle.nb_functions() << '\n';
        	std::cout << "Shape function d'indice 0 : " << shape_f_triangle.evaluate(0,point_ref) << '\n';
        	std::cout << "Shape function d'indice 1 : " << shape_f_triangle.evaluate(1,point_ref) << '\n';
        	std::cout << "Shape function d'indice 2 : " << shape_f_triangle.evaluate(2,point_ref) << '\n';
        	std::cout << "Vecteur gradient de la shape function d'indice 0 : " << shape_f_triangle.evaluate_grad(0,point_ref).x << " " << shape_f_triangle.evaluate_grad(0,point_ref).y << '\n';
        	std::cout << "Vecteur gradient de la shape function d'indice 1 : " << shape_f_triangle.evaluate_grad(1,point_ref).x << " " << shape_f_triangle.evaluate_grad(1,point_ref).y << '\n';
        	std::cout << "Vecteur gradient de la shape function d'indice 2 : " << shape_f_triangle.evaluate_grad(2,point_ref).x << " " << shape_f_triangle.evaluate_grad(2,point_ref).y << '\n';
        }
        
        void test_Ke()
        {
        	Mesh mesh;
        	mesh.load("data/square.mesh"); //récupère le cas square
        	ElementMapping mapping(mesh, false, 4);
        	ShapeFunctions shape_f_triangle(2,1);
        	
        	// Fixer un nombre dans quadrature : test avec 2
        	std::cout << "Test de Ke et K avec quadrature de 2" << '\n';
        	Quadrature quad = Quadrature::get_quadrature(2, false);
        	// test pour Ke en local
        	DenseMatrix Ke;
        	assemble_elementary_matrix(mapping, shape_f_triangle, quad, Simu::unit_fct, Ke);
        	std::cout << "Ke : " << '\n';
        	Ke.print();
        	
        	// test pour K en global
        	SparseMatrix K = SparseMatrix(mesh.nb_vertices());
        	local_to_global_matrix(mesh, 4, Ke, K);
        	std::cout << "K : " << '\n';
        	K.print();
        	
        	// Fixer un nombre dans quadrature : test avec 0
        	std::cout << "Test de Ke et K avec quadrature de 0" << '\n';
        	Quadrature quad_zero = Quadrature::get_quadrature(0, false);
       
        	DenseMatrix Ke_zero;
        	assemble_elementary_matrix(mapping, shape_f_triangle, quad_zero, Simu::unit_fct, Ke_zero);
        	std::cout << "Ke_zero  : " << '\n';
        	Ke_zero.print();
        	
        	SparseMatrix K_zero = SparseMatrix(mesh.nb_vertices());
        	local_to_global_matrix(mesh, 4, Ke_zero, K_zero);
        	std::cout << "K_zero : " << '\n';
        	K_zero.print();
        }
        
        void test_src_term()
        {
        	Mesh mesh;
        	mesh.load("data/square.mesh"); //récupère le cas square
        	ElementMapping mapping(mesh, false, 4);
        	ShapeFunctions shape_f_triangle(2,1);
        	
        	// Fixer un nombre dans quadrature : test avec 2
        	std::cout << "Test de Fe et F avec quadrature de 2" << '\n';
        	Quadrature quad = Quadrature::get_quadrature(2, false);
        	
        	// test pour Fe en local
        	std::vector< double > Fe(shape_f_triangle.nb_functions(), 0.);
        	assemble_elementary_vector(mapping, shape_f_triangle, quad, Simu::unit_fct, Fe);
        	for( int i = 0; i < shape_f_triangle.nb_functions(); ++i) {
        		std::cout << Fe[i] << " " << '\n';
        	}
        	
        	// test pour F en global
        	std::vector< double > F(mesh.nb_vertices(), 0.);
        	local_to_global_vector(mesh, false, 4, Fe, F);
        	for( int i = 0; i < mesh.nb_vertices(); ++i) {
        		std::cout << i << " " <<F[i] << '\n';
        	}
        	
        	// Fixer un nombre dans quadrature : test avec 0
        	std::cout << "Test de Fe et F avec quadrature de 0" << '\n';
        	Quadrature quad_zero = Quadrature::get_quadrature(0, false);
        	
        	// test pour Fe en local
        	std::vector< double > Fe_zero(shape_f_triangle.nb_functions(), 0.);
        	assemble_elementary_vector(mapping, shape_f_triangle, quad_zero, Simu::unit_fct, Fe_zero);
        	for( int i = 0; i < shape_f_triangle.nb_functions(); ++i) {
        		std::cout << Fe_zero[i] << " " << '\n';
        	}
        	
        	// test pour F en global
        	std::vector< double > F_zero(mesh.nb_vertices(), 0.);
        	local_to_global_vector(mesh, false, 4, Fe_zero, F_zero);
        	for( int i = 0; i < mesh.nb_vertices(); ++i) {
        		std::cout << i << " " << F_zero[i] << '\n';
        	}
        }

    }
}
