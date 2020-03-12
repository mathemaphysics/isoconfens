#include <iostream>
#include <fstream>
#include <getopt.h>
#include <gromacs/topology/topology.h>
#include <gromacs/fileio/trrio.h>
#include <boost/program_options.hpp>

int main(int argc, char **argv)
{
    /* Variables that need to be set */
    std::string trrFileName;
    gmx_int64_t frameStep;
    real frameTime;
    real frameLambda;
    rvec frameBox;
    int frameNumAtoms;
    rvec framePosition, frameVelocity, frameForce;

    /* Parse command line options */
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce help message")
        ("trr", po::value<std::string>(), "Set the TRR input file name")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("trr"))
    {
        std::cout << "TRR trajectory file set to "
                  << vm["trr"].as<std::string>() << std::endl;
        trrFileName = vm["trr"].as<std::string>();
    }
    else
        std::cout << "No input trajectory file name was set" << std::endl;

    /* Load the TRR file */
    t_fileio *trrPointer = gmx_trr_open(trrFileName.c_str(), "r");
    if (trrPointer != NULL)
    {
        while (gmx_trr_read_frame(
                trrPointer, &frameStep, &frameTime,
                &frameLambda, &frameBox, &frameNumAtoms,
                &framePosition, &frameVelocity, &frameForce
            )
        )
        {

        }
    }

    return 0;
}