#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A
{

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i )
        {
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

    const double triangle_P0[3] =
    {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] =
    {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] =
    {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] =
    {
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

    const double segment_P0[2] =
    {
        1., 0.5
    };

    const double segment_P2[4] =
    {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border )
        {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        }
        else if ( order == 2 && !border )
        {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        }
        else if ( order == 4 && !border )
        {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        }
        else if ( order == 6 && !border )
        {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        }
        else if ( order == 0 && border )
        {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        }
        else if ( order == 2 && border )
        {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        }
        else
        {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i )
        {
            if ( !border )
            {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            }
            else
            {
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
            : border_( border )
    {
        std::cout << "Test du constructeur ElementMapping pour l'element" << i << "  "<< std::endl;
        std::cout << "Coordonnees des sommets :"<<std::endl;
        if ( border ) std::cout << "(border)";
        // TODO
        // On détermine si l'élément est une frontière ou un triangle
        // on extrait les sommets de l'élément i du maillage M
        if (border)
        {
            // Supposons que les bords sont constitués de 2 sommets
            for (int v = 0; v<2; ++v) vertices_.push_back(M.get_edge_vertex(i, v));
            for (int v = 0; v<2; ++v) std::cout << vertices_[v].x<<"  "<<vertices_[v].y<<std::endl;
        }
        else
        {
            // Les triangles sont constitués de 3 sommets
            for (int v = 0; v<3; ++v) vertices_.push_back(M.get_triangle_vertex(i, v));
            for (int v = 0; v<3; ++v) std::cout << vertices_[v].x<<"  "<<vertices_[v].y<<std::endl;
        }
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';
        vertex r;
        // TODO
        if (border_)
        {
            double xi = x_r.x;
            r.x = (1-xi)*vertices_[0].x + xi*vertices_[1].x;
            r.y = (1-xi)*vertices_[0].y + xi*vertices_[1].y;
        }
        else
        {
            double xi = x_r.x ;
            double eta = x_r.y ;
            r.x = vertices_[0].x * (1 - xi - eta) + vertices_[1].x * xi + vertices_[2].x * eta;
            r.y = vertices_[0].y * (1 - xi - eta) + vertices_[1].y * xi + vertices_[2].y * eta;
        }
        /*
        // Supposons une interpolation affine pour un triangle
        if (!border_) {
           result.x = vertices_[0].x * x_r.x + vertices_[1].x * x_r.y + vertices_[2].x * (1 - x_r.x - x_r.y);
           result.y = vertices_[0].y * x_r.x + vertices_[1].y * x_r.y + vertices_[2].y * (1 - x_r.x - x_r.y);
        }*/
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        DenseMatrix J ;
        if (border_)
        {
            double dx_dxi = vertices_[1].x - vertices_[0].x;
            double dy_dxi = -vertices_[0].y + vertices_[1].y;
            J.set(0, 0, dx_dxi);
            J.set(1, 0, dy_dxi);
        }
        else
        {
            J.set_size(2,2);
            // Calcul des dérivées partielles pour la transformation d'un triangle
            double dx_dxi = vertices_[1].x - vertices_[0].x;   // ∂x/∂ξ = xB - xA
            double dx_deta = vertices_[2].x - vertices_[0].x;  // ∂x/∂η = xC - xA
            double dy_dxi = vertices_[1].y - vertices_[0].y;   // ∂y/∂ξ = yB - yA
            double dy_deta = vertices_[2].y - vertices_[0].y;  // ∂y/∂η = yC - yA
            J.set(0, 0, dx_dxi);
            J.set(0, 1, dx_deta);
            J.set(1, 0, dy_dxi);
            J.set(1, 1, dy_deta);
        }
        return J;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        // TODO
        DenseMatrix J = jacobian_matrix(x_r);
        // Calculer et retourner le déterminant de J
        if (border_)
        {
            double det_J = J.get(0,0)*J.get(0,0)+J.get(1,0)*J.get(1,0);
            return std::sqrt(det_J);
        }
        else
        {
            return J.det_2x2();
        }


        //return J.det_2x2();
        //return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
            : dim_( dim ), order_( order )
    {
        bool SF_construct = true;
        if (dim != 1 && dim != 2)
        {
            std::cout << "ShapeFunctions are only implemented in 1D or 2D" << std::endl;
            SF_construct = false;
        }
        if (order !=1)
        {
            std::cout << "Only order -1 ShapeFunctions are implemented" << std::endl;
            SF_construct = false ;
        }
        assert(SF_construct);
    }

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        // TODO
        //int nombre_fonction ;
        if (dim_ == 1)
        {
            return 2;
        }
        return 3;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        // TODO
        if (dim_ == 1)
        {
            double xi = x_r.x;
            switch (i)
            {
            case(0):
                            return 1 - xi;
            case (1):
                            return xi;
            }
        }
        else
{
            double xi = x_r.x;
            double eta = x_r.y;
            switch (i)
            {
            case (0): return 1 - xi - eta;
            case (1): return xi;
            case (2): return eta;
            }
        }
        return 0. ; // should not be reached

    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        // TODO

        //on peut faire avec switch mais on doit mettre break
        vec2 g ;
        if (dim_ == 1)
        {
            if (i == 0)
            {
                g.x = -1;
            }
            if (i == 1)
            {
                g.x = 1;
                g.y = 0;
            }
        }
        if (dim_ == 2)
        {
            if (i == 0)
            {
                g.x = -1;
                g.y = -1;
            }
            if (i == 1)
            {
                g.x = 1;
                g.y = 0;
            }
            if (i == 2)
            {
                g.x = 0;
                g.y = 1;
            }
        }

        /*
        double dx_dxi = vertices_[1].x - vertices_[0].x;   // ∂x/∂ξ = xB - xA
        double dx_deta = vertices_[2].x - vertices_[0].x;  // ∂x/∂η = xC - xA
        double dy_dxi = vertices_[1].y - vertices_[0].y;   // ∂y/∂ξ = yB - yA
        double dy_deta = vertices_[2].y - vertices_[0].y;  // ∂y/∂η = yC - yA
        */
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
        Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());

        for (int i = 0; i < reference_functions.nb_functions(); ++i)
        {
            for (int j = 0; j < reference_functions.nb_functions(); ++j)
            {
                Ke.set(i,j,0.);
                for (int q =0; q < quadrature.nb_points(); ++q)
                {
                    vertex p_q = quadrature.point(q);
                    double weight_q = quadrature.weight(q);
                    double detJ = elt_mapping.jacobian(p_q);
                    double k_value = coefficient(elt_mapping.transform(p_q));
                    DenseMatrix inv_J = elt_mapping.jacobian_matrix(p_q).invert_2x2();
                    vec2 grad_phi_i = inv_J.transpose().mult_2x2_2(reference_functions.evaluate_grad(i, p_q));
                    vec2 grad_phi_j = inv_J.transpose().mult_2x2_2(reference_functions.evaluate_grad(j, p_q));
                    double prod_scal = dot (grad_phi_i, grad_phi_j);
                    double integral_value = prod_scal * detJ * weight_q * k_value;
                    Ke.add(i,j,integral_value);
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
        //K=size
        //boucle sur les Ke
        //for (int k=0; k<M.nb_triangles(); ++k) {
        for (int local_i = 0; local_i<Ke.height(); ++local_i)
        {
            int global_i = M.get_triangle_vertex_index( t, local_i );
            for (int local_j=0; local_j<Ke.width(); ++local_j)
            {
                int global_j = M.get_triangle_vertex_index( t, local_j );
                double value = Ke.get(local_i,local_j );
                K.add(global_i, global_j,value);
            }
        }
        //}
    }


    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
        Fe.resize(reference_functions.nb_functions(), 0.0);

        for (int i = 0; i < reference_functions.nb_functions(); ++i)
    {
        //Fe.set(i,0.);
            for (int q =0; q < quadrature.nb_points(); ++q)
            {
                vertex p_q = quadrature.point(q);
                double weight_q = quadrature.weight(q);
                double detJ = elt_mapping.jacobian(p_q);

                double f_value = source(elt_mapping.transform(p_q));

                double phi_i = reference_functions.evaluate(i, p_q);

                double integral_value = weight_q * phi_i * f_value * detJ;
                Fe[i] += integral_value;
            }
        }
    }

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
        // TODO
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';

        //double large_penalty = std::numeric_limits<double>::max() / 1e8; // Une grande valeur mais pas trop pour éviter l'overflow

        // Supposer que M fournit un moyen d'itérer sur les attributs des bords et les nœuds associés
        /*for (int i = 0; i < M.nb_vertices(); ++i) {
            if (attribute_is_dirichlet[i]) {
                // Pour chaque nœud avec condition de Dirichlet
                for (int j = 0; j < K.nb_columns(); ++j) {
                    K.set(i, j, 0.0);  // Mise à zéro de la ligne entière
                    K.set(j, i, 0.0);  // Mise à zéro de la colonne entière
                }
                K.set(i, i, large_penalty);  // Mettre une grande valeur sur la diagonale
                F[i] = values[i] * large_penalty;  // Ajuster le vecteur F
            }
        }
        */

        std::vector< bool > processed_vertices(values.size(), false);
        double penalty_coefficient = 10000.;
        for ( int edge = 0; edge < M.nb_edges(); edge++ )
        {
            int edge_attribute = M.get_edge_attribute(edge);
            if ( attribute_is_dirichlet[edge_attribute] )
            {
                for ( int v = 0; v < 2; v++ )
                {
                    int vertex_index = M.get_edge_vertex_index( edge, v);
                    if ( !processed_vertices[vertex_index] )
                    {
                        processed_vertices[vertex_index] = true;
                        K.add(vertex_index, vertex_index, penalty_coefficient);
                        F[vertex_index] += penalty_coefficient*values[vertex_index];
                    }
                }
            }
        }
        std::cout << "Dirichlet conditions applied.\n";
    }

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
