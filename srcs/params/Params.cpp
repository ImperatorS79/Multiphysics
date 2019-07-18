/**
 * \file Params.cpp
 * \brief Implementation of the required function to load a SolverParams struct from file.
 */

#include <iostream>
#include <fstream>
#include <cerrno> //Change with something more c++ ;-)
#include <cstring>
#include <vector>
#include "nlohmann/json.hpp"
#include "Params.hpp"
#include "../physics/boundaryConditions.hpp"
#include "../physics/fluxes.hpp"
#include "../physics/phiPsis.hpp"
#include "../physics/sources.hpp"


/**
 * \brief Load boundary conditions from a file
 * \param j JSON object describing the parameters.
 * \param solverParams The structure in which the parameters are loaded.
 * \param fileName The name of the parameters file which has been opened (for debug output).
 * \return true if the loading succeeds, false otherwise.
 */
static bool handleBoundaryCondition(const nlohmann::json& j,
                                    SolverParams& solverParams,
                                    const std::string& fileName)
{
    auto initBoundaryConditions = j["physics"]["initialBoundaryConditions"];
    unsigned int nBC = 0;
    for(auto initBoundaryCondition : initBoundaryConditions)
    {
        ibc tempCondition;
        bool error = false;
        if(solverParams.problemType == "transport")
        {
            if(initBoundaryCondition["type"] == "constant")
                tempCondition.ibcFunc = constant;

            else if(initBoundaryCondition["type"] == "sinusTransport")
                tempCondition.ibcFunc = sinusTransport;

            else if(initBoundaryCondition["type"] == "gaussianTransport")
                tempCondition.ibcFunc = gaussianTransport;

            else if(initBoundaryCondition["type"] == "freeTransport")
                tempCondition.ibcFunc = freeTransport;

            else if(initBoundaryCondition["type"] == "gaussian2DTransport")
                tempCondition.ibcFunc = gaussian2DTransport;

            else
                error = true;
        }
        else if(solverParams.problemType == "shallow")
        {
            if(initBoundaryCondition["type"] == "constant")
                tempCondition.ibcFunc = constant;

            else if(initBoundaryCondition["type"] == "affineShallow")
                tempCondition.ibcFunc = affineShallow;

            else if(initBoundaryCondition["type"] == "sinusShallow")
                tempCondition.ibcFunc = sinusShallow;

            else if(initBoundaryCondition["type"] == "sinusAffShallow")
                tempCondition.ibcFunc = sinusAffShallow;

            else if(initBoundaryCondition["type"] == "reflectShallow")
                tempCondition.ibcFunc = reflectShallow;

            else if(initBoundaryCondition["type"] == "gaussian2DShallow")
                tempCondition.ibcFunc = gaussian2DShallow;

            else if(initBoundaryCondition["type"] == "gaussian1DShallowX")
                tempCondition.ibcFunc = gaussian1DShallowX;

            else if(initBoundaryCondition["type"] == "gaussian1DShallowY")
                tempCondition.ibcFunc = gaussian1DShallowY;

            else if(initBoundaryCondition["type"] == "openShallow")
                tempCondition.ibcFunc = openShallow;

            else if(initBoundaryCondition["type"] == "openAffShallow")
                tempCondition.ibcFunc = openAffShallow;

            else
                error = true;

        }
        else if(solverParams.problemType == "shallowLin")
        {
            if(initBoundaryCondition["type"] == "constant")
                tempCondition.ibcFunc = constant;

            else if(initBoundaryCondition["type"] == "sinusShallowLin")
                tempCondition.ibcFunc = sinusShallowLin;

            else if(initBoundaryCondition["type"] == "reflectShallowLin")
                tempCondition.ibcFunc = reflectShallowLin;

            else if(initBoundaryCondition["type"] == "gaussian2DShallowLin")
                tempCondition.ibcFunc = gaussian2DShallowLin;

            else if(initBoundaryCondition["type"] == "gaussian1DShallowXLin")
                tempCondition.ibcFunc = gaussian1DShallowXLin;

            else if(initBoundaryCondition["type"] == "gaussian1DShallowYLin")
                tempCondition.ibcFunc = gaussian1DShallowYLin;

            else if(initBoundaryCondition["type"] == "openShallowLin")
                tempCondition.ibcFunc = openShallowLin;

            else
                error = true;
        }
        else if(solverParams.problemType == "acousticLin")
        {
            if(initBoundaryCondition["type"] == "constant")
                tempCondition.ibcFunc = constant;

            else if(initBoundaryCondition["type"] == "sinusAcousticLin")
                tempCondition.ibcFunc = sinusAcousticLin;

            else if(initBoundaryCondition["type"] == "reflectAcousticLin")
                tempCondition.ibcFunc = reflectAcousticLin;

            else if(initBoundaryCondition["type"] == "gaussian2DAcousticLin")
                tempCondition.ibcFunc = gaussian2DAcousticLin;

            else if(initBoundaryCondition["type"] == "gaussian1DAcousticLinX")
                tempCondition.ibcFunc = gaussian1DAcousticLinX;

            else if(initBoundaryCondition["type"] == "gaussian1DAcousticLinY")
                tempCondition.ibcFunc = gaussian1DAcousticLinY;

            else if(initBoundaryCondition["type"] == "openAcousticLin")
                tempCondition.ibcFunc = openAcousticLin;

            else
                error = true;
        }
        if(error)
        {
            std::cerr << "Unhandled boundary condition type "
                      << initBoundaryCondition["type"]
                      << " for boundary " << std::endl
                      << initBoundaryCondition["physicalGroup"] << " for problem type "
                      << solverParams.problemType << " in parameter file "
                      << fileName <<std::endl;

            return false;
        }

        tempCondition.coefficients = initBoundaryCondition["coefficients"].get<std::vector<double>>();

        if(initBoundaryCondition["physicalGroup"] == "Init_Cond")
            solverParams.initCondition = tempCondition;

        else
        {
            solverParams.boundaryConditions[
                                initBoundaryCondition["physicalGroup"]] = tempCondition;

            ++nBC;
        }
    }

    std::cout << "Initial condition present and " << nBC
              << " boundary conditions present in "<< std::endl
              << "file " << fileName << std::endl;

    return true;
}

/**
 * \brief Load general parameters from a file
 * \param j JSON object describing the parameters.
 * \param solverParams The structure in which the parameters are loaded.
 * \param fileName The name of the parameters file which has been opened (for debug output).
 * \return true if the loading succeeds, false otherwise.
 */
static bool loadGeneralParams(const nlohmann::json& j,
                              const std::string& fileName,
                              SolverParams& solverParams)
{
    std::string temp = j["general"]["spaceIntegrationType"];
    if(temp.compare(0, 5, "Gauss") != 0
       || !(temp.substr(5, temp.size() - 5).find_first_not_of("0123456789")
            == std::string::npos))
    {
        std::cerr << "Unexpected space integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        return false;
    }
    solverParams.spaceIntType = temp;

    temp = j["general"]["basisFunctionType"];
    if(!(temp == "Lagrange" || temp == "Isoparametric"))
    {
        std::cerr << "Unexpected basis function type " << temp
                  << " in parameter file " << fileName << std::endl;

        return false;
    }
    solverParams.basisFuncType = temp;

    temp = j["general"]["timeIntegrationType"];
    if(!(temp == "RK1" || temp == "RK2" || temp == "RK3" || temp == "RK4"))
    {
        std::cerr << "Unexpected time integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        return false;
    }
    solverParams.timeIntType = temp;

    temp = j["general"]["solverType"];
    if(!(temp == "strong" || temp == "weak"))
    {
        std::cerr << "Unexpected solver type " << temp
                  << " in parameter file " << fileName << std::endl;

        return false;
    }
    solverParams.solverType = temp;

    solverParams.simTime = j["general"]["simulationTime"];

    solverParams.timeStep = j["general"]["simulationTimeSteps"];

    solverParams.simTimeDtWrite = j["general"]["simulationTimeToWrite"];

    return true;
}

/**
 * \brief Load physical parameters from a file
 * \param j JSON object describing the parameters.
 * \param solverParams The structure in which the parameters are loaded.
 * \param fileName The name of the parameters file which has been opened (for debug output).
 * \return true if the loading succeeds, false otherwise.
 */
static bool loadPhysicsParams(const nlohmann::json& j,
                              const std::string& fileName,
                              SolverParams& solverParams)
{
    std::string temp;

    solverParams.problemType = j["physics"]["problemType"];
    if(solverParams.problemType == "shallow")
    {
        solverParams.nUnknowns = 3;
        solverParams.flux = fluxShallow;
    }
    else if(solverParams.problemType == "transport")
    {
        solverParams.nUnknowns = 1;
        solverParams.flux = fluxTransport;
    }
    else if(solverParams.problemType == "shallowLin")
    {
        solverParams.nUnknowns = 3;
        solverParams.flux = fluxShallowLin;
    }
    else if(solverParams.problemType == "AcousticLin")
    {
        solverParams.nUnknowns = 3;
        solverParams.flux = fluxAcousticLin;
    }
    else
    {
        std::cerr << "Unexpected problem type " << temp
                  << " in parameter file " << fileName << std::endl;

        return false;
    }

    bool error = false;
    std::vector<std::string> whatToWrite = j["physics"]["whatToWrite"];
    if(solverParams.problemType == "shallow" ||
               solverParams.problemType == "shallowLin")
    {
        solverParams.whatToWrite.resize(5);
        solverParams.viewTags.resize(5);
        std::fill(solverParams.whatToWrite.begin(),
                  solverParams.whatToWrite.end(), false);
        for (unsigned short i = 0 ; i < whatToWrite.size() ; ++i)
        {
            if(whatToWrite[i] == "H")
                solverParams.whatToWrite[0] = true;

            else if(whatToWrite[i] == "u")
                solverParams.whatToWrite[1] = true;

            else if(whatToWrite[i] == "v")
                solverParams.whatToWrite[2] = true;

            else if(whatToWrite[i] == "sKE")
                solverParams.whatToWrite[3] = true;

            else if(whatToWrite[i] == "vField")
                solverParams.whatToWrite[4] = true;

            else
                error = true;
        }
    }
    else if(solverParams.problemType == "AcousticLin")
    {
        solverParams.whatToWrite.resize(5);
        solverParams.viewTags.resize(5);
        std::fill(solverParams.whatToWrite.begin(),
                  solverParams.whatToWrite.end(), false);
        for (unsigned short i = 0 ; i < whatToWrite.size() ; ++i)
        {
            if(whatToWrite[i] == "p'")
                solverParams.whatToWrite[0] = true;

            else if(whatToWrite[i] == "u'")
                solverParams.whatToWrite[1] = true;

            else if(whatToWrite[i] == "v'")
                solverParams.whatToWrite[2] = true;

            else if(whatToWrite[i] == "sKE'")
                solverParams.whatToWrite[3] = true;

            else if(whatToWrite[i] == "vField'")
                solverParams.whatToWrite[4] = true;

            else
                error = true;
        }
    }
    else if(solverParams.problemType == "transport")
    {
        solverParams.whatToWrite.resize(1);
        solverParams.viewTags.resize(1);
        std::fill(solverParams.whatToWrite.begin(),
                  solverParams.whatToWrite.end(), false);
        for (unsigned short i = 0 ; i < whatToWrite.size() ; ++i)
        {
            if(whatToWrite[i] == "u")
                solverParams.whatToWrite[0] = true;

            else
                error = true;
        }
    }
    if(error)
    {
        std::cerr << "Unexpected writing parameters ("
                  << j["physics"]["whatToWrite"].dump() << ") for problem type "
                  << solverParams.problemType << std::endl;

        return false;
    }

    if(solverParams.problemType == "shallow")
        solverParams.write = writeShallow;

    else if(solverParams.problemType == "shallowLin")
        solverParams.write = writeShallowLin;

    else if(solverParams.problemType == "transport")
        solverParams.write = writeTransport;

    else if(solverParams.problemType == "AcousticLin")
        solverParams.write = writeAcousticLin;

    temp = j["physics"]["numericalFlux"];
    if(solverParams.problemType == "shallow")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFShallow;
        }
        else if(temp == "Roe")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = Roe;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "shallowLin")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFShallowLin;
        }
        else if(temp == "Roe")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = RoeLin;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "AcousticLin")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFAcousticLin;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "transport")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFTransport;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    if(error)
    {
        std::cerr << "Unexpected flux type type " << temp
                  << " for problem type " << solverParams.problemType
                  << " in parameter file " << fileName << std::endl;

        return false;
    }

    error = false;
    solverParams.fluxCoeffs = j["physics"]["fluxCoefficients"].get<std::vector<double>>();
    if(solverParams.problemType == "shallow")
    {
        if(solverParams.fluxCoeffs.size() !=1)
            error = true;
    }
    else if(solverParams.problemType == "shallowLin")
    {
        if(solverParams.fluxCoeffs.size() != 2)
            error = true;

    }
    else if(solverParams.problemType == "transport")
    {
        if(solverParams.fluxCoeffs.size() != 2)
            error = true;

    }
    else if(solverParams.problemType == "AcousticLin")
    {
        if(solverParams.fluxCoeffs.size() != 4)
            error = true;

    }
    if(error)
    {
        std::cerr << "Unexpected number of flux coefficients ("
                  << solverParams.fluxCoeffs.size() << ") for problem type "
                  << solverParams.problemType << std::endl;

        return false;
    }

    temp = j["physics"]["sourceTerms"];
    solverParams.IsSourceTerms = true;
    error = false;
    if(solverParams.problemType == "shallow")
    {
        if(temp == "no")
        {
            solverParams.IsSourceTerms = false;
            solverParams.sourceType = temp;
        }
        else if(temp == "sourceShallowCstGradCstFrict"
                && solverParams.problemType == "shallow")
        {
            solverParams.sourceType = temp;
            solverParams.sourceTerm = sourceShallowCstGradCstFrict;
        }
        else if(temp == "sourceShallowCstGradQuadFrict"
                && solverParams.problemType == "shallow")
        {
            solverParams.sourceType = temp;
            solverParams.sourceTerm = sourceShallowCstGradQuadFrict;
        }
        else if(temp == "shallowLinCst" && solverParams.problemType == "shallowLin")
        {
            solverParams.sourceType = temp;
            solverParams.sourceTerm = sourceShallowLinCst;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "shallowLin" ||
            solverParams.problemType == "transport" ||
            solverParams.problemType == "AcousticLin")
    {
        solverParams.IsSourceTerms = false;
        solverParams.sourceType = "no";
    }

    if(error)
    {
            std::cerr << "Unexpected source function ("
                      << solverParams.sourceType << ") for problem type "
                      << solverParams.problemType << std::endl;

            return false;
    }


    error = false;
    solverParams.sourceCoeffs = j["physics"]["sourceCoefficients"].get<std::vector<double>>();
    if(solverParams.problemType == "shallow")
    {
        if(solverParams.sourceType == "sourceShallowCstGradCstFrict")
        {
            if(solverParams.IsSourceTerms)
            {
                if(solverParams.sourceCoeffs.size() != 5)
                    error = true;
            }
        }
        else if(solverParams.sourceType == "sourceShallowCstGradQuadFrict")
        {
            if(solverParams.IsSourceTerms)
            {
                if(solverParams.sourceCoeffs.size() != 5)
                    error = true;
            }
        }
        else if(solverParams.sourceType == "shallowLinCst")
        {
            if(solverParams.IsSourceTerms)
            {
                if(solverParams.sourceCoeffs.size() != 1)
                    error = true;
            }
        }
    }

    if(error)
    {
        std::cerr << "Unexpected number of source terms coefficients ("
                      << solverParams.sourceCoeffs.size() << ") for problem type "
                      << solverParams.problemType << std::endl;

        return false;
    }

    if(!handleBoundaryCondition(j, solverParams, fileName))
        return false;

    return true;
}

//Documentation in .hpp
bool loadSolverParams(const std::string& fileName, SolverParams& solverParams)
{
    std::ifstream paramFile(fileName);

    if(!paramFile.is_open())
    {
        std::cerr << "Something went wrong when trying to read the file "
                  << fileName << std::endl
                  << std::strerror(errno) << std::endl;

        paramFile.close();
        return false;
    }

    nlohmann::json j;
    paramFile >> j;
    paramFile.close();

    if(!loadGeneralParams(j, fileName, solverParams))
        return false;

    if(!loadPhysicsParams(j, fileName, solverParams))
        return false;

    // display the parameters
    std::cout   << "Number of Gauss points: " << solverParams.spaceIntType
                << std::endl
                << "Type of basis function: " << solverParams.basisFuncType
                << std::endl
                << "Time intgeration scheme: " << solverParams.timeIntType
                << std::endl
                << "Formulation type: " << solverParams.solverType
                << std::endl
                << "Simulation time duration: " << solverParams.simTime << "s"
                << std::endl
                << "Time step: " << solverParams.timeStep << "s"
                << std::endl
                << "Time between data writing: " << solverParams.simTimeDtWrite
                << "s"
                << std::endl
                << "Problem type: " << solverParams.problemType
                << std::endl
                << "Source terms: " << solverParams.sourceType
                << std::endl
                << "Numerical Flux: " << solverParams.fluxType
                << std::endl;

    return true;
}
