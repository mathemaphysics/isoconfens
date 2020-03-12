#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <gromacs/topology/topology.h>
#include <gromacs/fileio/trrio.h>
#include <boost/program_options.hpp>

int main(int argc, char **argv)
{
    /* Variables that need to be set */
    std::string trrFileName;
    int frameNumAtoms;
    gmx_int64_t frameStep;
    real frameTime, frameLambda;
    rvec frameBox, framePosition, frameVelocity, frameForce;

    /* Parse command line options */
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce help message")
        ("trr", po::value<std::string>(), "Set the TRR input file name")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    /* Notify the user of variables set via command line options */
    po::notify(vm);

    /* If help is given then show the help and exit */
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }

    /* If TRR file is given notify user that was heard */
    if (vm.count("trr"))
    {
        std::cout << "TRR trajectory file set to "
                  << vm["trr"].as<std::string>() << std::endl;
        trrFileName = vm["trr"].as<std::string>();
    }
    else
    {
        std::cout << "No input trajectory file name was set" << std::endl;
        return 1;
    }

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
            /* Periodically notify which step */
            //if (frameStep % 100 == 0)
            //    std::cout << "--> Read step " << std::setw(7) << frameStep << std::endl;
        }
    }

    return 0;
}