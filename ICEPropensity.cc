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
    int frameNumAtoms, bufferLevel;
    u_int64_t frameStep;
    real frameTime, frameLambda, *msdAccumulator;
    rvec frameBox, *frameBuffer;
    gmx_trr_header_t frameHeader;
    gmx_bool outOk;
    const int frameBufferSize = 100;
    const int maxTimeDiff = 100;

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

    bool firstTrajectory = true;
    while (!trrFiles.empty())
    {
        std::string file;
        file = trrFiles.top();
        std::cout << "Running trajectory " << file << std::endl;

        /* Load the TRR file */
        t_fileio *trrPointer = gmx_trr_open(file.c_str(), "r");
        if (trrPointer != NULL)
        {
            if (gmx_trr_read_frame_header(trrPointer, &frameHeader, &outOk))
            {
                std::cout << "Header has been read" << std::endl;
                std::cout << "There are " << std::setw(7) << frameHeader.natoms << " in the system" << std::endl;

                /* Allocate only the first time */
                if (firstTrajectory)
                {
                    frameBuffer = new rvec[frameBufferSize*frameHeader.natoms];
                    msdAccumulator = new real[frameHeader.natoms];
                    firstTrajectory = false;
                }
            }
            else
                throw(std::exception());

            bufferLevel = 0;
            frameStep = 0;
            while (
                gmx_trr_read_frame_data(
                    trrPointer, &frameHeader, &frameBox,
                    &frameBuffer[bufferLevel * frameHeader.natoms],
                    NULL, NULL
                )
            )
            {
                if (bufferLevel >= frameBufferSize - 1)
                {
                    int origin = frameStep;
                    for (int i = 0; i < frameBufferSize; i++)
                    {
                        for (int j = 0; j < maxTimeDiff; j++)
                        {
                            for (int k = 0; k < frameHeader.natoms; k++)
                            {
                                rvec diffVector;
                                rvec_sub(frameBuffer[i * frameHeader.natoms + k],
                                         frameBuffer[(i + j) * frameHeader.natoms + k],
                                         diffVector);
                                real metricDistance = norm2(diffVector);
                                msdAccumulator[k] = msdAccumulator[k] + metricDistance;
                            }
                        }
                    }
                    bufferLevel = 0;
                }
                else
                    bufferLevel = bufferLevel + 1;
                
                frameStep = frameStep + 1;
            }
        }
        gmx_trr_close(trrPointer);
        trrFiles.pop();
    }

    /* Only deallocate if you allocated, i.e. found at least a trajectory */
    if (!firstTrajectory)
    {
        delete[] frameBuffer;
        delete[] msdAccumulator;
    }

    return 0;
}
