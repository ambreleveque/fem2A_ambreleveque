#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border ) //constructeur d'ElementMapping
    {
    	if (border) std::cout << "(border)"; // s'il y a border alors segment
    	std::cout << '\n';
    	
    	if (border) { // cas d'un segment donc max que deux vertices = deux points
    		for (int v_local_index = 0; v_local_index < 2; v_local_index++) { //on note v le vertex local index
    			vertices_.push_back(M.get_edge_vertex(i,v_local_index));
    		}
    	}
    	else { //cas d'un triangle donc trois vertices
    		for (int v_local_index = 0; v_local_index < 3; v_local_index++) {
    			vertices_.push_back(M.get_triangle_vertex(i,v_local_index));
    		}
    	}
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        std::cout << "[ElementMapping] transform reference to world space " << '\n';
        
        vertex r ; // dans le réel
        if (border_) { //cas segment
        	r.x = (1 - x_r.x) * vertices_[0].x + x_r.x * vertices_[1].x;
        	r.y = (1 - x_r.x) * vertices_[0].y + x_r.x * vertices_[1].y;
        	std::cout << "Coordonnées du vertice du segment dans le réel " << r.x << " " << r.y << '\n';
        }
        else { //cas triangle
        	r.x = (1 - x_r.x - x_r.y)* vertices_[0].x + x_r.x * vertices_[1].x + x_r.y * vertices_[2].x;
        	r.y = (1 - x_r.x - x_r.y)* vertices_[0].y + x_r.x * vertices_[1].y + x_r.y * vertices_[2].y;
        	std::cout << "Coordonnées du vertice du triangle dans le réel " << r.x << " " << r.y << '\n';
        }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian matrix " << '\n';
 
        DenseMatrix J ;
        if (border_) {
        	J.set_size(2,1);
        	J.set(0,0, -vertices_[0].x + vertices_[1].x);
        	J.set(1,0, -vertices_[0].y + vertices_[1].y);
        }
        else {
        	J.set_size(2,2);
        	J.set(0,0,vertices_[1].x - vertices_[0].x);
        	J.set(1,0,vertices_[1].y - vertices_[0].y);
        	J.set(0,1,vertices_[2].x - vertices_[0].x);
        	J.set(1,1,vertices_[2].y - vertices_[0].y);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian determinant " << '\n';
 
        DenseMatrix J = jacobian_matrix(x_r);
        if (border_) {
        	DenseMatrix T = J.transpose();
        	float produit = J.get(0,0)*T.get(0,0) + J.get(1,0)*T.get(0,1);
        	float det = sqrt(produit);
        	std::cout << "Le determinant est : " << det << '\n';
        	return det;
        }
        else {
        	float det = J.det_2x2();
        	std::cout << "Le determinant est : " << det << '\n';
        	return det;
        }
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        if (dim_ != 1 && dim != 2) {
        	std::cout << "Attention, vous avez entré une mauvaise dimension" << '\n';
        }
        if (order_ != 1) {
        	std::cout << "Attention, vous avez entré un ordre supérieur à 1" << '\n';
        }
    }

    int ShapeFunctions::nb_functions() const
    {
        std::cout << "[ShapeFunctions] number of functions" << '\n';
        if (dim_ == 1) {
        	return 1;
	}
        else {
        	return 3;
	}
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        if (dim_ ==1) {
        	switch(i) {
        		case 0 :
        			return 1 - x_r.x ;
        		break;
        		case 1 :
        			return x_r.x ;
        		break;
        	}
        }		
        else {
        	switch(i) {
        		case 0 :
        			return 1 - x_r.x - x_r.y ; break;
        		case 1 :
        			return x_r.x ; break;
        		case 2 :
        			return x_r.y ; break;
        	}
        }
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        vec2 g ;
        
        if (dim_==1) {
        	switch(i) {
        		case 0 :
        			g.x = -1 ; break;
        		case 1 :
        			g.x = 1 ; break;
        	}
        }
        else {
        	switch(i) {
        		case 0 :
        			g.x = -1 ; g.y = -1 ; break;
        		case 1 :
        			g.x = 1 ; g.y = 0 ; break;
        		case 2 :
        			g.x = 0 ; g.y = 1 ; break;
        	}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        std::cout << "compute elementary matrix" << '\n';
        // taille Ke est le nbre de points d'interpolation
        Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());
        for (int i=0; i < reference_functions.nb_functions(); ++i) {
        	for (int j = 0; j < reference_functions.nb_functions(); ++j) {
        		Ke.set(i, j, 0.);
        		for (int k = 0; k < quadrature.nb_points(); ++k) {
        			// points de gauss = points de la quadrature sur lesquels la somme de Ke est réalisée
        			vertex ptg_q = quadrature.point(k);
        			DenseMatrix J = elt_mapping.jacobian_matrix(ptg_q);
        			DenseMatrix inv_J = J.invert_2x2().transpose();
        			vec2 grad_i = reference_functions.evaluate_grad(i, ptg_q);
        			vec2 grad_j = reference_functions.evaluate_grad(j, ptg_q);
        			Ke.add(i, j, quadrature.weight(k) * coefficient(elt_mapping.transform(ptg_q)) * dot(inv_J.mult_2x2_2(grad_i), inv_J.mult_2x2_2(grad_j)) * elt_mapping.jacobian(ptg_q));
        		}
        	}
        }
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        std::cout << "Ke -> K" << '\n';
        // taille de K est le nbre de points d'interpolation globale, ie nbre de points du maillage
        for (int ligne = 0; ligne < Ke.height(); ++ligne){
        	// parcours de la matrice Ke sur ses lignes et colonnes et récupération des indices
        	int i = M.get_triangle_vertex_index(t, ligne);
        	for(int colonne = 0; colonne < Ke.width(); ++colonne){
        		int j = M.get_triangle_vertex_index(t, colonne);
        		K.add(i, j, Ke.get(ligne, colonne));
        	}
        }
    }
    
    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        for (int i = 0; i < reference_functions.nb_functions(); ++i) {
        	for (int k = 0; k < quadrature.nb_points(); ++k) {
        		const vertex ptg_q = quadrature.point(k);
        		//On redimensionne Fe et on la rempli avec 0 au début puis avec les valeurs de la boucle
        		Fe[i] += quadrature.weight(k) * source(elt_mapping.transform(ptg_q)) * reference_functions.evaluate(i, ptg_q) * elt_mapping.jacobian(ptg_q);
        		std::cout << "Calcul de la somme" << '\n';
        	}
        }
    }
    
// Condition de Neumann non realisee
    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        
        // condition d'un segment
        if (border) {
        	// parcours des lignes de la matrice Fe pour récupérer les index
        	for ( int ligne = 0; ligne < Fe.size(); ++ligne) {
        		F[M.get_edge_vertex_index(i, ligne)] += Fe[ligne];
        	}
        }
        // condition d'un triangle
        else {
        	// parcours des lignes de la matrice Fe
        	for ( int ligne = 0; ligne < Fe.size(); ++ligne){
        		F[M.get_triangle_vertex_index(i, ligne)] += Fe[ligne];
        	}
        }
    }
    
    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        std::vector<bool> processed_vertices(values.size(), false);
        double penalty_coefficient = 10000.;
        
       // parcours de chaque segment
        for (int edge = 0; edge < M.nb_edges(); ++edge) {
        	int edge_attribute = M.get_edge_attribute(edge);
        	if ( attribute_is_dirichlet[edge_attribute] ) {
        		// parcours de chaque noued de chaque segment
        		for (int n = 0; n < 2; ++n){
      				// application de la condition de Dirichlet
        			int vertex_index = M.get_edge_vertex_index(edge,n);
        			// condition de vérification qu'un noeud est traité qu'une seule fois
        			if ( !processed_vertices[vertex_index] ) {
        				processed_vertices[vertex_index] = true;
        				K.add(vertex_index, vertex_index, penalty_coefficient);
        				F[vertex_index] += penalty_coefficient * values[vertex_index];
        			}
        		}
        	}
        }
    }

// Non realisee
    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
