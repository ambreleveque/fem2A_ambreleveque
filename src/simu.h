#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            Mesh mesh;
            mesh.load(mesh_filename);
            SparseMatrix K(mesh.nb_vertices());
            std::vector< double > F(mesh.nb_vertices(), 0.);
            
            // parcours des triangles consituant le maillage
            for ( int triangle = 0; triangle < mesh.nb_triangles(); ++triangle) {
            	ElementMapping mapping(mesh, false, triangle);
            	ShapeFunctions shape_f_triangle(2,1);
            	Quadrature quad = Quadrature::get_quadrature(2);
            	DenseMatrix Ke;
            	// on utilise unit_fct pour le calcul de Ke car k = 1
            	assemble_elementary_matrix(mapping, shape_f_triangle, quad, unit_fct, Ke);
            	local_to_global_matrix(mesh, triangle, Ke, K);
            }
            
            // condition de Dirichlet
            std::vector< double > values(mesh.nb_vertices());
            std::vector< bool > attribut_dirichlet(2, false);
            attribut_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            
            for (int i = 0; i < mesh.nb_vertices(); ++i) {
            	values[i] = xy_fct(mesh.get_vertex(i));
            }
            apply_dirichlet_boundary_conditions(mesh, attribut_dirichlet,values, K, F);
            
            // résolution du système linéaire
            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            // sauvergarde
            std::string export_name ="pure_dirichlet";
            mesh.save(export_name+".mesh"); /* sauvergarde du maillage */
            save_solution(u, export_name+".bb"); /* sauvergarde de la solution du pb */
        }
	
	void dirichlet_with_src_pb(const std::string& mesh_filename, bool verbose)
	{
            std::cout << "Solving a Dirichlet problem with a source term" << std::endl;
            Mesh mesh;
            mesh.load(mesh_filename);
            SparseMatrix K(mesh.nb_vertices());
            std::vector< double > F(mesh.nb_vertices(), 0.);
            
            // parcours des triangles consituant le maillage
            for ( int triangle = 0; triangle < mesh.nb_triangles(); ++triangle) {
            	ElementMapping mapping(mesh, false, triangle);
            	ShapeFunctions shape_f_triangle(2,1);
            	Quadrature quad = Quadrature::get_quadrature(2);
            	// on utilise unit_fct pour le calcul de Ke car k = 1
            	DenseMatrix Ke;
            	assemble_elementary_matrix(mapping, shape_f_triangle, quad, unit_fct, Ke);
            	local_to_global_matrix(mesh, triangle, Ke, K);
            	std::vector< double > Fe(shape_f_triangle.nb_functions(), 0.);
            	assemble_elementary_vector(mapping, shape_f_triangle, quad, unit_fct, Fe);
            	local_to_global_vector(mesh, false, triangle, Fe, F);
            }
            // Condition de Dirichlet
            std::vector< double > values(mesh.nb_vertices());
            std::vector< bool > attribut_dirichlet(2, false);
            attribut_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            
            for (int i = 0; i < mesh.nb_vertices(); ++i) {
            	values[i] = zero_fct(mesh.get_vertex(i));
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribut_dirichlet, values, K, F);
            
            // solve du système linéaire
            std::vector< double > u(mesh.nb_vertices());
            solve(K, F, u);
            
            // sauvegarde
            std::string export_name = "dirichlet_with_source_term";
            mesh.save(export_name+".mesh"); /* sauvergarde du maillage */
            save_solution(u, export_name+".bb"); /* sauvergarde de la solution du pb */
	}
    }

}
