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

        bool test_quadrature(int ordre, bool bord)
        {
        	std::cout << "test Quadrature" << std::endl ;
        	Quadrature Q;
        	Q = Quadrature::get_quadrature(ordre, bord);
        	std::cout << Q.nb_points()<< std::endl;
        	double sum=0;
        	for (int i=0; i<Q.nb_points(); ++i) {
        		std::cout << Q.point(i).x << " " << Q.point(i).y << std::endl;
        		std::cout << Q.weight(i) << std::endl;
        		sum = sum + Q.weight(i);
        		}
        		std::cout << sum << std::endl;
 		return true;
    	}


    	bool test_constructeur_elementmapping (int elementIndex, bool Border) {
    	    Mesh mesh;
            mesh.load("data/square.mesh");

            // on peut distiguer les triangles et edge ? à améliorer
            // Création de l'objet ElementMapping
            ElementMapping mapping(mesh, Border, elementIndex);
            std::cout << "Mapping du point (xi=0.2, eta=0.4) : " << '\n';
            vertex test_ver = {0.2,0.4};
            test_ver = mapping.transform(test_ver);
            std::cout << test_ver.x<<"  "<<test_ver.y<< '\n';

            DenseMatrix test_J;

            test_J = mapping.jacobian_matrix(test_ver);
            std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
            test_J.print();

            double det_J = mapping.jacobian(test_ver);
            std::cout << det_J << std::endl;

            //double test_det =

    	    return true;
    	    }

        bool test_shapefunction () {
    	    Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};
            std::cout << "Test du constructeur ShapeFunctions\n";
            // Test avec une dimension de 1 et un ordre correct de 1
            ShapeFunctions Segments(1, 1);
            int nombre_fonction = Segments.nb_functions();
            std::cout<<nombre_fonction<<std::endl;
            assert(Segments.nb_functions() == 2);  // Un segment devrait avoir 2 fonctions de forme
            std::cout << "Passed: Line segment with correct order.\n";

            // Test avec une dimension de 2 et un ordre correct de 1
            ShapeFunctions Triangle(2, 1);
            assert(Triangle.nb_functions() == 3);  // Un triangle devrait avoir 3 fonctions de forme
            std::cout << "Passed: Triangle with correct order.\n";

            double resultat_eval = Triangle.evaluate(1,test_ver);
            std::cout << resultat_eval<<std::endl;

            vec2 vecteur = Triangle.evaluate_grad( 2, test_ver );
            std::cout << vecteur.x <<"  "<<vecteur.y<<std::endl;
    	    return true;
    	    }

            double unit_fct( vertex v ) {
                return 1.;
            }




    	    bool test_assemble_elementary_matrix (int elementIndex, bool Border) {
            std::cout << "Testing assemble_elementary_matrix\n";

    	    Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};

              // Création de l'objet ElementMapping
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            // Création de l'objet shape function d'un triangle linéaire
            ShapeFunctions reference_functions( 2, 1 );
            // Création de l'objet quadrature
            Quadrature Q;
        	Q = Quadrature::get_quadrature(2, false);

            //Création de la dense matrix
            DenseMatrix Ke;
            Ke.set_size(3,3);

            //fonction coefficient
            double coefficient_k = unit_fct( test_ver );

            //test de la fonction

            std::cout<<"tout va bien"<<coefficient_k<<std::endl;

            assemble_elementary_matrix(elt_mapping,reference_functions,Q,unit_fct,Ke);

            Ke.print();
    	    return true;
    	    }



    	    bool test_local_to_global_matrix (int elementIndex, bool Border) {
            std::cout << "Testing local_to_global_matrix\n";

    	    Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};

              // Création de l'objet ElementMapping
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            // Création de l'objet shape function d'un triangle linéaire
            ShapeFunctions reference_functions( 2, 1 );
            // Création de l'objet quadrature
            Quadrature Q;
        	Q = Quadrature::get_quadrature(2, false);

        	double coefficient_k = unit_fct( test_ver );

            //Création de la dense matrix
            DenseMatrix Ke;
            Ke.set_size(3,3);
            assemble_elementary_matrix(elt_mapping,reference_functions,Q,unit_fct,Ke);

            Ke.print();
            SparseMatrix K(mesh.nb_triangles());


            //test de la fonction
            int t=4;

            local_to_global_matrix(mesh, t,Ke, K);
            K.print();
            std::cout<<"tout va bien"<<coefficient_k<<std::endl;
    	    return true;
    	    }

            bool test_assemble_elementary_vector (int elementIndex, bool Border){
            std::cout << "Testing assemble_elementary_matrix\n";

    	    Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};

              // Création de l'objet ElementMapping
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            // Création de l'objet shape function d'un triangle linéaire
            ShapeFunctions reference_functions( 2, 1 );
            // Création de l'objet quadrature
            Quadrature Q;
        	Q = Quadrature::get_quadrature(2, false);

            //Création du vecteur Fe
            std::vector< double > Fe;

            //fonction source
            double coefficient_f = unit_fct( test_ver );

            //test de la fonction

            std::cout<<"tout va bien fe ?"<<coefficient_f<<std::endl;

            assemble_elementary_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);

            for (int i = 0; i < Fe.size(); ++i) {
            std::cout << Fe[i] <<std::endl;
            }
    	    return true;
    	    }


    	    bool test_local_to_global_vector (int elementIndex, bool Border) {
            std::cout << "Testing local_to_global_vector\n";

    	    Mesh mesh;
            mesh.load("data/square.mesh");

            vertex test_ver = {0.2,0.4};

              // Création de l'objet ElementMapping
            ElementMapping elt_mapping(mesh, Border, elementIndex);
            // Création de l'objet shape function d'un triangle linéaire
            ShapeFunctions reference_functions( 2, 1 );
            // Création de l'objet quadrature
            Quadrature Q;
        	Q = Quadrature::get_quadrature(2, false);

            //Création de la matrice Fe
            std::vector< double > Fe;
            assemble_elementary_vector(elt_mapping,reference_functions,Q,unit_fct,Fe);
            std::vector< double > F;
            F.resize(mesh.nb_triangles(),0.0);
            int t = elementIndex;
            local_to_global_vector(mesh, Border, t,Fe, F);

            std::cout<<"K ? ?"<<std::endl;

            for (int i = 0; i < F.size(); ++i) {
            std::cout << F[i] <<std::endl;
            }
    	    return true;
    	    }

    	}
}
