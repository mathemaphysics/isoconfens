#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <gromacs/topology/topology.h>
#include <gromacs/fileio/trrio.h>
#include <gromacs/math/vectypes.h>
#include <gromacs/math/vec.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>

int main(int argc, char **argv)
{
    /* Variables that need to be set */
    std::string trrFileBaseName;
    int frameNumAtoms;
    u_int64_t frameStep;
    real frameTime, frameLambda;
    rvec frameBox, tempVector, *firstFramePosition, *framePosition, *frameVelocity, *frameForce;
    gmx_trr_header_t frameHeader;
    gmx_bool outOk;

    /* Set up namespaces for convenience */
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;
    namespace bs = boost;

    /* Parse command line options */
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce help message")
        ("trr", po::value<std::string>(), "Set the TRR input file name")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
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
        std::cout << "TRR trajectory file base name set to "
                  << vm["trr"].as<std::string>() << std::endl;
        trrFileBaseName = vm["trr"].as<std::string>();
    }
    else
    {
        std::cout << "No input trajectory file name was set" << std::endl;
        return 1;
    }

    /* Set up iteration through files in the current directory */
    std::stack<std::string> trrFiles;
    fs::path curdir = fs::path(".");
    fs::directory_iterator endItr;
    bs::regex fnameRegex{trrFileBaseName + "[0-9]{3}\\.trr"};
    for (fs::directory_iterator itr(curdir); itr != endItr; ++itr)
    {
        if (fs::is_regular_file(itr->path()))
        {
            bs::cmatch match;
            std::string curTraj = itr->path().string();
            if (bs::regex_search(curTraj.c_str(), match, fnameRegex))
                trrFiles.push(curTraj);
        }
    }

    while (!trrFiles.empty())
    {
        std::string file;
        file = trrFiles.top();
        std::cout << "Running trajectory " << file << std::endl;
        trrFiles.pop();
    }

    return 0;

    /* Load the TRR file */
    t_fileio *trrPointer = gmx_trr_open(trrFileBaseName.c_str(), "r");
    if (trrPointer != NULL)
    {
        if (gmx_trr_read_frame_header(trrPointer, &frameHeader, &outOk))
        {
            std::cout << "Header has been read" << std::endl;
            std::cout << "There are " << std::setw(7) << frameHeader.natoms << " in the system" << std::endl;
            std::cout << "The system size: " << std::setw(7) << frameHeader.step << std::endl;

            firstFramePosition = new rvec[frameHeader.natoms];
            framePosition = new rvec[frameHeader.natoms];
            frameVelocity = new rvec[frameHeader.natoms];
            frameForce = new rvec[frameHeader.natoms];
        }
        frameStep = 0;

        while (trrPointer)
        {
            gmx_trr_read_frame_data(trrPointer, &frameHeader, &frameBox,
                                       framePosition, frameVelocity, frameForce);
            int origin = frameStep;
            //for (int i = 0; i < frameHeader.natoms; i++)
            //{
            //    rvec_sub(framePosition[i], firstFramePosition[i], tempVector);
            //    norm2(tempVector);
            //}
            frameStep = frameStep + 1;
        }
    }

    return 0;
}
