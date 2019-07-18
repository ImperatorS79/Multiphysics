#include <iostream>
#include "flux.hpp"
#include <iostream>

// see .hpp file for description
void fluxAcousticLin(Field& field, PartialField& partialField,
                      const SolverParams& solverParams, bool boundary)
{
    // the physical flux is given by
    //  Fx = [rho*co^2*u' + u0*p',     p'/rho + u'*u0,            u'*v0]
    //  Fy = [rho*co^2*v' + v0*p',              v'*u0,   p'/rho + v'*v0]

    if(boundary)
    {
        // flux for the (x, y) coordinates for p'
        partialField.FluxAtBC[0][0] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                                    * solverParams.fluxCoeffs[1]*partialField.uAtBC[1]
                                    + solverParams.fluxCoeffs[2]*partialField.uAtBC[0];
        partialField.FluxAtBC[1][0] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                                    * solverParams.fluxCoeffs[1]*partialField.uAtBC[2]
                                    + solverParams.fluxCoeffs[3]*partialField.uAtBC[0];

        // flux for the (x, y) coordinates for u'
        partialField.FluxAtBC[0][1] = partialField.uAtBC[0]/solverParams.fluxCoeffs[0]
                                    + solverParams.fluxCoeffs[2]*partialField.uAtBC[1];
        partialField.FluxAtBC[0][2] = solverParams.fluxCoeffs[3]*partialField.uAtBC[1];

        // flux for the (x, y) coordinates for v'
        partialField.FluxAtBC[1][1] = solverParams.fluxCoeffs[2]*partialField.uAtBC[2];
        partialField.FluxAtBC[1][2] = partialField.uAtBC[0]/solverParams.fluxCoeffs[0]
                                    + solverParams.fluxCoeffs[3]*partialField.uAtBC[2];
    }
    else
    {
        // flux for the (x, y) coordinates for H
        field.flux[0][0] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                         * solverParams.fluxCoeffs[1]*field.u[1]
                         + solverParams.fluxCoeffs[2]*field.u[0];
        field.flux[1][0] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                         * solverParams.fluxCoeffs[1]*field.u[2]
                         + solverParams.fluxCoeffs[3]*field.u[0];

        // flux for the (x, y) coordinates for u
        field.flux[0][1] = field.u[0]/solverParams.fluxCoeffs[0]
                         + solverParams.fluxCoeffs[2]*field.u[1];
        field.flux[0][2] = solverParams.fluxCoeffs[3]*field.u[1];

        // flux for the (x, y) coordinates for v
        field.flux[1][1] = solverParams.fluxCoeffs[2]*field.u[2];
        field.flux[1][2] = field.u[0]/solverParams.fluxCoeffs[0]
                         + solverParams.fluxCoeffs[3]*field.u[2];
    }
}
