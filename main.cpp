#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_test_quadrature = false;
    const bool t_test_element_mapping = false;
    const bool t_test_shapefunctions = false;
    const bool t_test_assemblary_matrix = false;
    const bool t_test_local_to_global_matrix = false;
    const bool t_test_apply_dirichlet_boundary_conditions = false;
    const bool t_test_assemblary_vector = false;
    const bool t_test_local_to_global_vector = false;
    const bool t_test_assemblary_neumann_vector = true;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_test_quadrature) Tests::test_quadrature(2,false); //Test of the quadrature of a triangle to order 2
    if( t_test_element_mapping) Tests::test_elementmapping(4, true);//Test de la classe Element Mapping pour un segment
    if( t_test_element_mapping) Tests::test_elementmapping(4, false);//Test de la classe Element Mapping pour un triangle
    if (t_test_shapefunctions) Tests::test_shapefunctions(0);// Test of the 0-th shape function
    if (t_test_assemblary_matrix) Tests::test_assemble_elementary_matrix (4, false);//Test of assemble_elementary_matrix
    if (t_test_local_to_global_matrix) Tests::test_local_to_global_matrix (4, false);
    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );
    if (t_test_apply_dirichlet_boundary_conditions) Tests::test_apply_dirichlet_boundary_conditions ("data/square.mesh", verbose);
    if (t_test_assemblary_vector) Tests::test_assemble_elementary_vector (4, false);
    if (t_test_local_to_global_vector) Tests::test_local_to_global_vector (4, false);
    if (t_test_local_to_global_vector) Tests::test_local_to_global_vector (4, true);
    if (t_test_assemblary_neumann_vector) Tests::test_assemble_elementary_neumann_vector (4, true);

}

void run_simu()
{

    const bool s_pure_dirichlet_pb = false;
    const bool s_source_dirichlet_pb = false;
    const bool s_sinus_bump_dirichlet_pb = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( s_pure_dirichlet_pb ) {
        Simu::pure_dirichlet_pb("data/square_fine.mesh", verbose);
    }
    if( s_source_dirichlet_pb ) {
        Simu::source_dirichlet_pb("data/square.mesh", verbose);
    }
    if( s_sinus_bump_dirichlet_pb ) {
        Simu::sinus_bump_dirichlet_pb("data/square_fine.mesh", verbose);
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
