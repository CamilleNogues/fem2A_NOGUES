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

namespace FEM2A
{
    namespace Tests
    {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for ( int v = 0; v < mesh.nb_vertices(); v++ )
            {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for ( int ed = 0; ed < mesh.nb_edges(); ed++ )
            {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                << mesh.get_edge_vertex_index(ed, 1) << " "
                << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for ( int tr = 0; tr < mesh.nb_triangles(); tr++ )
            {
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
        /**START OF IMPLEMENTATION TESTS**/

        bool test_quadrature(int ordre, bool Border)
        {
            std::cout << "[TEST OF THE QUADRATURE]" << std::endl ;
            Quadrature Q;
            Q = Quadrature::get_quadrature(ordre, Border);
            if (Border == false)
            {
                std::cout << "Quadrature d'un triangle" << std::endl;
                std::cout << "Nombre de points :" << Q.nb_points()<< std::endl;
            }
            else
            {
                std::cout << "Quadrature d'un segment" <<  std::endl;
                std::cout << "Nombre de points :" << Q.nb_points()<< std::endl;
            }
            double sum = 0;
            for (int i=0; i<Q.nb_points(); ++i)
            {
                std::cout << "coordonnees X du point " << i << " : " << Q.point(i).x << std::endl;
                std::cout << "coordonnees Y du point " << i << " : " << Q.point(i).y << std::endl;
                std::cout << " poids du point " << i << " : " << Q.weight(i) << std::endl;
                sum = sum + Q.weight(i);
            }
            std::cout << "Valeur de la quadrature a l'ordre " << ordre << " : " << sum << std::endl;
            return true;
        }


        bool test_elementmapping (int elementIndex, bool Border)
        {
            if (Border == true)
            {
                std::cout << "[ELEMENT MAPPING TEST FOR AN EDGE]"<< std::endl;
            }
            else
            {
                std::cout << "[ELEMENT MAPPING TEST FOR A TRIANGLE]"<< std::endl;
            }
            Mesh mesh;
            mesh.load("data/square.mesh");
            /**
            * \brief Creation of the ElementMapping object
            */
            ElementMapping mapping(mesh, Border, elementIndex); //Test of the constructor directly in the src\fem.cpp
            /**
            * \brief Test of the transform function which transforms a point from local space to global space
            */
            std::cout << "- Transform function test - " << std::endl;
            std::cout << "Mapping of the point (xi=0.2, eta=0.4) : " << std::endl;
            vertex point_test = {0.2,0.4};
            point_test = mapping.transform(point_test);
            std::cout << point_test.x<<"  "<<point_test.y<< std::endl;
            /**
            * \brief Test of the jacobian matrix function which computes the jacobian matrix of the mapping.
            */
            DenseMatrix test_J;
            test_J = mapping.jacobian_matrix(point_test);
            std::cout << "Jacobian matrix at the point (xi=0.2, eta=0.4) : " << std::endl;
            test_J.print();
            double det_J = mapping.jacobian(point_test);
            std::cout << "determinant = " <<det_J << std::endl;

            return true;
        }

        bool test_shapefunctions (int shapefunctionIndex)
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            vertex point_test = {0.2,0.4};
            int i = shapefunctionIndex; //indice of the shape function
            std::cout << "[SHAPE FUNCTIONS TEST]"<< std::endl;
            /*Construction des elements segment et triangle à l'ordre 1*/
            ShapeFunctions Segment(1, 1);
            ShapeFunctions Triangle(2, 1);
            /**Test de la méthode nb_functions.
            *Un segment devrait avoir 2 fonctions de forme
            *Un triangle devrait avoir 3 fonctions de forme
            **/
            int nombre_fonction_segment = Segment.nb_functions();
            int nombre_fonction_triangle = Triangle.nb_functions();
            std::cout<<"Number of shape functions for the edge :" <<nombre_fonction_segment<<std::endl;
            std::cout<<"Number of shape functions for the triangle :" <<nombre_fonction_triangle<<std::endl;

            assert(Segment.nb_functions() == 2);  // Passed: Line segment with correct order.
            assert(Triangle.nb_functions() == 3);  // Passed: Triangle with correct order.

            double res_evaluate_segment = Segment.evaluate(i,point_test);
            double res_evaluate_triangle = Triangle.evaluate(i,point_test);

            std::cout << "Evaluate result of an edge for the " << i << "-th shape function : " << res_evaluate_segment<<std::endl;
            std::cout << "Evaluate result of a triangle for the " << i << "-th shape function : "<< res_evaluate_triangle<<std::endl;

            vec2 vecteur = Triangle.evaluate_grad( i, point_test);
            std::cout << "Evaluate gradient result of a triangle for the " << i << "-th shape function : " << vecteur.x <<"  "<<vecteur.y<<std::endl;

            return true;
        }


        double unit_fct( vertex v )
        {
            return 1.;
        }




        bool test_assemble_elementary_matrix (int elementIndex, bool Border)
        {
            std::cout << "[ASSEMBLE ELEMENTARY MATRIX TEST]" << std::endl;
            Mesh mesh;
            mesh.load("data/square.mesh");
            /**
            * \brief Creation of the ElementMapping object
            */
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            /**
            * \brief Creation of the shape function object of a linear triangle
            */
            ShapeFunctions reference_functions( 2, 1 );
            /**
            * \brief Creation of the quadrature object
            */
            Quadrature Q;
            Q = Quadrature::get_quadrature(2, false);
            /**
            * \brief Creation of the dense matrix Ke
            */
            DenseMatrix Ke;
            Ke.set_size(3,3);
            /**
            * \brief Test of the assemble_elementary_matrix method which computes the elementary matrix Ke associated
            * to a triangle defined by its ElementMapping
            */
            assemble_elementary_matrix(elt_mapping,reference_functions,Q,unit_fct,Ke);
            std::cout << "Elementary matrix Ke : " << '\n';
            Ke.print();
            return true;
        }


        bool test_local_to_global_matrix (int elementIndex, bool Border)
        {
            std::cout << "[LOCAL TO GLOBAL MATRIX TEST]" << std::endl;
            Mesh mesh;
            mesh.load("data/square.mesh");
            /**
            * \brief Creation of the ElementMapping object
            */
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            /**
            * \brief Creation of the shape function object of a linear triangle
            */
            ShapeFunctions reference_functions( 2, 1 );
            /**
            * \brief Creation of the quadrature object
            */
            Quadrature Q;
            Q = Quadrature::get_quadrature(2, false);
            /**
            * \brief Creation of the dense matrix Ke
            */
            DenseMatrix Ke;
            Ke.set_size(3,3);
            assemble_elementary_matrix(elt_mapping,reference_functions,Q,unit_fct,Ke);
            Ke.print();
            /**
            * \brief Creation of the sparse matrix K
            */
            SparseMatrix K(mesh.nb_triangles());
            /**
            * \brief Test of the local_to_global_matrix method which adds the contribution Ke of triangle t to
            *        the global matrix K.
            */
            local_to_global_matrix(mesh, elementIndex,Ke, K);
            std::cout << "Ke -> K" << '\n';
            for (int local_j=0; local_j<3; ++local_j)
            {
                int global_j = mesh.get_triangle_vertex_index( elementIndex, local_j );
                double value = Ke.get(local_j,local_j );
                K.add(global_j, global_j,value);
                std::cout<<"Le numero global du point "<<local_j<<" est "<<global_j<<std::endl;
            }
            K.print();
            return true;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        bool test_apply_dirichlet_boundary_conditions (const std::string& mesh_filename, bool verbose)
        {
            std::cout << "[APPLY DIRICHLET BOUNDARY CONDITIONS TEST]" << std::endl;
            Mesh M;
            M.load(mesh_filename);

            /**
            * \brief Creation of the sparse matrix K
            */
            SparseMatrix K_glob(M.nb_vertices());
            std::vector<double> F_glob(M.nb_vertices(),0.);
            /**
            * \brief Test of the local_to_global_matrix method which adds the contribution Ke of triangle t to
            *        the global matrix K.
            */
            std::cout << "Solving a pure Dirichlet problem" << std::endl;

            for (int t=0; t<M.nb_triangles(); ++t)
            {
                ElementMapping my_map(M, false, t);
                ShapeFunctions my_shpfct(2,1);
                Quadrature my_quad = Quadrature::get_quadrature(2);
                DenseMatrix Ke;
                assemble_elementary_matrix(my_map, my_shpfct, my_quad, unit_fct, Ke);
                local_to_global_matrix(M,t,Ke,K_glob);
            }


            std::vector<bool> att_is_dirichlet(2,false);
            att_is_dirichlet[1] = true;
            M.set_attribute(unit_fct, 1, true);
            std::vector<double> imposed_values(M.nb_vertices());
            for (int i=0; i<M.nb_vertices(); ++i)
            {
                imposed_values[i] = xy_fct(M.get_vertex(i));
            }
            apply_dirichlet_boundary_conditions(M, att_is_dirichlet, imposed_values, K_glob, F_glob);
            K_glob.print();
            return true;
        }




        bool test_assemble_elementary_vector (int elementIndex, bool Border)
        {
            std::cout << "[ASSEMBLE ELEMENTARY VECTOR TEST]" << std::endl;
            Mesh mesh;
            mesh.load("data/square.mesh");
            /**
            * \brief Creation of the ElementMapping object
            */
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            /**
            * \brief Creation of the shape function object of a linear triangle
            */
            ShapeFunctions reference_functions( 2, 1 );
            /**
            * \brief Creation of the quadrature object
            */
            Quadrature Q;
            Q = Quadrature::get_quadrature(2, false);
            /**
            * \brief Creation of the vector Fe
            */
            std::vector< double > Fe;
            /**
            * \brief Test of the assemble_elementary_vector method which computes the elementary vector Fe
            *        associated to a triangle defined by its ElementMapping due to the source term.
            */
            assemble_elementary_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);
            std::cout << "Elementary vector Fe : " << '\n';
            for (int i = 0; i < Fe.size(); ++i)
            {
                std::cout << Fe[i] <<std::endl;
            }
            return true;
        }


        bool test_local_to_global_vector (int elementIndex, bool Border)
        {
            std::cout << "[LOCAL TO GLOBAL VECTOR TEST]" << std::endl;
            Mesh mesh;
            mesh.load("data/square.mesh");
            /**
            * \brief Creation of the ElementMapping object
            */
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            /**
            * \brief Creation of the shape function object of a linear triangle (dim = 2)
            *        or an edge (dim = 1).
            */
            int dim; //
            if ( Border == true)
            {
                dim = 1;
            }
            else
            {
                dim = 2;
            }
            ShapeFunctions reference_functions( dim, 1 );
            /**
            * \brief Creation of the quadrature object
            */
            Quadrature Q;
            Q = Quadrature::get_quadrature(2, Border);
            /**
            * \brief Creation of the vector Fe
            */
            std::vector< double > Fe;
            assemble_elementary_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);
            /**
            * \brief Creation of the matrix F
            */
            std::vector< double > F;
            if (Border)
            {
                assemble_elementary_neumann_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);
                F.resize(mesh.nb_edges(),0.0);
            }
            else
            {
                assemble_elementary_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);
                F.resize(mesh.nb_triangles(),0.0);
            }
            local_to_global_vector(mesh, Border, elementIndex,Fe, F);

            return true;
        }

        bool test_assemble_elementary_neumann_vector (int elementIndex, bool Border)
        {
            std::cout << "[ASSEMBLE ELEMENTARY NEUMANN VECTOR TEST]" <<std::endl;
            Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};

            // Création de l'objet ElementMapping
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            // Création de l'objet shape function d'un segment linéaire
            ShapeFunctions reference_functions( 1, 1 );
            // Création de l'objet quadrature
            Quadrature Q;
            Q = Quadrature::get_quadrature(2, Border);

            //Création du vecteur Fe
            std::vector< double > Fe;

            //fonction source
            double coefficient_f = unit_fct( test_ver );

            //test de la fonction

            std::cout<<"tout va bien fe ?"<<coefficient_f<<std::endl;
            std::cout << "compute elementary vector Fe (neumann condition)" << '\n';


            assemble_elementary_neumann_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);

            for (int i = 0; i < Fe.size(); ++i)
            {
                std::cout << Fe[i] <<std::endl;
            }
            return true;
        }
    }
}
