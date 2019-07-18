#include "phiPsi.hpp"

// see .hpp file for description
void LFAcousticLin(const Edge& edge, Field& field, PartialField& partialField,
                unsigned int j, double factor, bool boundary, unsigned int indexJ,
                unsigned int indexFrontJ, const SolverParams& solverParams)
{
    // speed of sound parameter
    double c0 = solverParams.fluxCoeffs[1];

    double lambdaIn =   (field.u[1][indexJ]*edge.normal[0]
                    +   field.u[2][indexJ]*edge.normal[1]);

    lambdaIn = (lambdaIn >= 0) ? lambdaIn + c0 : -lambdaIn + c0;

    double lambdaOut;

    // compute the numerical flux
    if(boundary)
    {
        // computation of the absolute value of the eigenvalue outside the element
        lambdaOut = (partialField.uAtBC[1]*edge.normal[0]
                    + partialField.uAtBC[2]*edge.normal[1]);

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + c0 : -lambdaOut + c0;

        // computation of the value of C in the LF scheme
        double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

        //computation of g
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                        + partialField.FluxAtBC[dim][unk]
                    + C*edge.normal[dim]*(field.u[unk][indexJ]
                        - partialField.uAtBC[unk]))/2;
            }
        }
    }
    else
    {
        // computation of the absolute value of the eigenvalue outside the element
        lambdaOut = (field.u[1][indexFrontJ]*edge.normal[0]
                  + field.u[2][indexFrontJ]*edge.normal[1]);

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + c0 : -lambdaOut + c0;

        // computation of the value of C in the LF scheme
        double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

        //computation of g
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                        + field.flux[dim][unk][indexFrontJ]
                    + C*edge.normal[dim]*(field.u[unk][indexJ]
                        - field.u[unk][indexFrontJ]))/2;
            }
        }
    }
}
