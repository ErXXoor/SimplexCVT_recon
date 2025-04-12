//
// Created by hongbo on 09/04/25.
//
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include "Base/MeshAdaptor.h"
#include "RVD/SurfRec.h"
int main(int argc, char **argv){
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("poly");
    GEO::CmdLine::import_arg_group("co3ne");

    std::string input_filename;
    std::string output_filename;
    int nb_pts = 1000;
    int dim = 6;

    if(argc == 5){
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        nb_pts = std::stoi(argv[4]);
    }else {
        std::cout << "Parameter length error" << std::endl;
        return 1;
    }

    GEO::Mesh point_in(dim);

    Base::MeshAdaptor::LoadHDXYZ(input_filename, point_in, dim);

    double R = GEO::bbox_diagonal(point_in);
//    GEO::mesh_repair(point_in, GEO::MESH_REPAIR_COLOCATE, 1e-6*R);

    double radius =5*0.01*R;

    RVD::SurfRec reconstructor;
    reconstructor.init(point_in);
    reconstructor.Co3ne_rec(output_filename);


    return 0;
}