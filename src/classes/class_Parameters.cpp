/*! \file class_parameter.cpp
 * \ingroup Parser
 * ### Created on 5/25/18 by Matthew Varga
 * ### Purpose
 * ***
 * Class member functions for class Parameter
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "classes/class_Membrane.hpp"
#include "classes/class_Parameters.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

double Parameters::dt = 0.0;
std::vector<long long int> Parameters::lastUpdateTransition {};

namespace {

long long int parse_integer_timestep_value(const std::string& value)
{
    std::size_t parsedLength { 0 };
    double parsedValue { 0.0 };

    try {
        parsedValue = std::stod(value, &parsedLength);
    } catch (const std::exception&) {
        throw std::invalid_argument("Not a valid timestep value.");
    }

    if (parsedLength != value.size() || !std::isfinite(parsedValue)
        || parsedValue < static_cast<double>(std::numeric_limits<long long int>::min())
        || parsedValue > static_cast<double>(std::numeric_limits<long long int>::max())
        || parsedValue != std::trunc(parsedValue)) {
        throw std::invalid_argument("Not a valid timestep value.");
    }

    return static_cast<long long int>(parsedValue);
}

void set_timestep_value(long long int& target, const std::string& value, const char* label)
{
    target = parse_integer_timestep_value(value);
    std::cout << "Read in " << label << ": " << target << " timeSteps" << std::endl;
}

} // namespace

// this is so we can compare parsed keywords with the enumerations (we need them as strings)
std::map<const std::string, ParamKeyword> parmKeywords = {
    { "nummoltypes", ParamKeyword::numMolTypes }, { "numtotalspecies", ParamKeyword::numTotalSpecies },
    { "nitr", ParamKeyword::nItr },
    { "fromrestart", ParamKeyword::fromRestart }, { "timewrite", ParamKeyword::timeWrite },
    { "trajwrite", ParamKeyword::trajWrite }, { "timestep", ParamKeyword::timeStep },
    { "assocdissocwrite", ParamKeyword::assocDissocWrite},
    { "numtotalcomplex", ParamKeyword::numTotalComplex },
    { "mass", ParamKeyword::mass }, { "restartwrite", ParamKeyword::restartWrite },
    { "pdbwrite", ParamKeyword::pdbWrite },
    { "bondedcomplexwrite", ParamKeyword::bondedComplexWrite },
    { "overlapseplimit", ParamKeyword::overlapSepLimit }, { "name", ParamKeyword::name },
    { "checkpoint", ParamKeyword::checkPoint }, { "scalemaxdisplace", ParamKeyword::scaleMaxDisplace },
    { "transitionwrite", ParamKeyword::transitionWrite }, 
    { "clusteroverlapcheck", ParamKeyword::clusterOverlapCheck }, 
    { "rngwrite", ParamKeyword::rngwrite },
};

void Parameters::set_value(std::string value, ParamKeyword keywords)
{
    /*! \ingroup Parser
     * \brief Sets the parameters based on the enumeration key.
     * @param value value of the parameter as a string
     * @param keywords the keyword parsed from the input file to match to the enumeration Keywords
     */
    try {
        auto key = static_cast<std::underlying_type<ParamKeyword>::type>(keywords);
        switch (key) {
        case 0:
            this->numMolTypes = std::stoi(value);
            std::cout << "Read in numMolTypes: " << this->numMolTypes << std::endl;
            break;
        case 1:
            this->numTotalSpecies = std::stoi(value);
            std::cout << "Read in numTotalSpecies: " << this->numTotalSpecies << std::endl;
            break;
        case 2:
            set_timestep_value(this->nItr, value, "nItr");
            break;
        case 3:
            this->fromRestart = read_boolean(value);
            std::cout << "Read in fromRestart: " << std::boolalpha << this->fromRestart << std::endl;
            break;
        case 4:
            set_timestep_value(this->timeWrite, value, "timeWrite");
            break;
        case 5:
            set_timestep_value(this->trajWrite, value, "trajWrite");
            break;
        case 6:
            this->timeStep = std::stod(value);
            std::cout << "Read in timeStep: " << this->timeStep << " us" << std::endl;
            break;
        case 7:
            this->numTotalComplex = std::stoi(value);
            std::cout << "Read in numTotalComplex: " << this->numTotalComplex << std::endl;
            break;
        case 8:
            this->mass = std::stod(value);
            std::cout << "Read in mass: " << this->mass << std::endl;
            break;
        case 10:
            set_timestep_value(this->restartWrite, value, "restartWrite");
            break;
        case 11:
            set_timestep_value(this->pdbWrite, value, "pdbWrite");
            break;
        case 12:
            this->overlapSepLimit = std::stod(value);
            std::cout << "Read in overlapSepLimit: " << this->overlapSepLimit << " nm" << std::endl;
            break;
        case 13:
            this->name = value;
            std::cout << "Read in name: " << value << std::endl;
            break;
        case 14:
            set_timestep_value(this->checkPoint, value, "checkPoint");
            break;
        case 15:
            this->scaleMaxDisplace = std::stod(value);
            std::cout << "Read in scaleMaxDisplace: " << this->scaleMaxDisplace << std::endl;
            break;
        case 16:
            set_timestep_value(this->transitionWrite, value, "transitionWrite");
            break;
        case 17:
            this->clusterOverlapCheck = read_boolean(value);
            std::cout << "Read in clusterOverlapCheck: " << std::boolalpha << this->clusterOverlapCheck << std::endl;
            break;
        case 18:
            this->assocDissocWrite = read_boolean(value);
            std::cout << "Read in assocDissocWrite: " << std::boolalpha << this->assocDissocWrite << std::endl;
            break;
        case 19:
            this->rngwrite = read_boolean(value);
            std::cout << "Read in RNGwrite: " << std::boolalpha << this->rngwrite << std::endl;
            break;
        case 20:
            set_timestep_value(this->bondedComplexWrite, value, "bondedComplexWrite");
            break;
        default:
            throw std::invalid_argument("Not a valid keyword.");
        }
    } catch (std::invalid_argument& e) {
        std::cout << e.what() << '\n';
        exit(1);
    }
}

void Parameters::parse_paramFile(std::ifstream& paramFile)
{
    /*! \ingroup Parser
     * \brief Main function to parse the parameters block of an input file
     */

    while (!paramFile.eof()) {
        std::string line;
        // see if we're at the end of the parameter block
        auto startPos = paramFile.tellg();
        getline(paramFile, line);

        // remove spaces
        line.erase(
            std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
        std::string tmpLine { line };
        std::transform(tmpLine.begin(), tmpLine.end(), tmpLine.begin(), ::tolower);

        if (tmpLine == "endparameters") {
            paramFile.seekg(startPos);
            return;
        }
        // if the line starts with a comment, skip it
        if (tmpLine[0] == '#')
            continue;
        else
            remove_comment(line);

        bool gotValue { false };
        std::string buffer;
        for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
            if (std::isalnum(*lineItr))
                buffer += std::tolower(static_cast<char>(*lineItr));
            else if (*lineItr == '=') {
                auto keyFind = parmKeywords.find(buffer);
                line.erase(line.begin(), lineItr + 1); // + 1 removes the '=' sign. could make this erase(remove_if)

                // find the value type from the keyword and then set that parameter
                if (keyFind != parmKeywords.end()) {
                    // std::cout << "Keyword found: " << keyFind->first << '\n';
                    this->set_value(line, keyFind->second);
                    gotValue = true;
                    break;
                } else {
                    std::cout << "Warning, ignoring unknown keyword " << buffer << '\n';
                    break;
                }
            }
        }
    }
}

void Parameters::display()
{
    std::cout << "Number of iterations: " << nItr << " timesteps\n";
    std::cout << "Timestep: " << timeStep << " us\n";
    std::cout << "Timestep log interval (timeWrite): " << timeWrite << " timesteps\n";
    std::cout << "Restart file write interval: " << restartWrite << " timesteps\n";
    std::cout << "Coordinate write interval (trajWrite): " << trajWrite << " timesteps\n";
    std::cout << "PDB Coordinate write interval: " << pdbWrite << " timesteps\n";
    std::cout << "Checkpoint write interval: " << checkPoint << " timesteps\n";
    std::cout << "overlapSepLimit: " << overlapSepLimit << " nm\n";
    std::cout << "Transition matrix write interval: " << transitionWrite << " timesteps\n";
    std::cout << "Bonded complex JSON write interval: " << bondedComplexWrite << " timesteps\n";
    std::cout << "ClusterOverlapCheck: " << clusterOverlapCheck << "\n";
    std::cout << "RNGwrite: " << rngwrite << "\n";

    std::cout << "Molecule specific parameters:\n";
    std::cout << "Number of unique molecule types: " << numMolTypes << '\n';
    std::cout << "Total number of unique interfaces and states, including product states: " << numTotalSpecies << '\n';
    std::cout << "Total number of complexes in system at start: " << numTotalComplex << '\n';
    std::cout << "Total number of units (molecules + interfaces) in system at start: " << numTotalUnits << '\n';
    std::cout << "Maximum allowed number of unique 2D reactions: " << max2DRxns << '\n';
}
