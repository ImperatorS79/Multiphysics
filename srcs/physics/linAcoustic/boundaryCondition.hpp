#ifndef linAcoustic_boundaryCondition_hpp_included
#define linAcoustic_boundaryCondition_hpp_included

#include <vector>
#include "../../solver/field.hpp"



/**
 * \brief Compute a wave of the shape A*sin(2*pi*nu*t + phi) + B -- for the
 * linear acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the sinus: coeffs[0] = A, coeffs[1] = nu,
 * coeffs[2] = phi, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void sinusAcousticLin(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a physical reflection -- for linear acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient (not used here).
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void reflectAcousticLin(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a physical opening -- for linear acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient (not used here).
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void openAcousticLin(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a 2D gaussian: A*exp(-(x-x0)^2/(2*var_y)-(y-y0)^2/(2*var_y)) -- for
 * linear acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = y_0, coeffs[4] = var_y.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian2DAcousticLin(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a 1D gaussian along x: A*exp(-(x-x0)^2/(2*var_x) + B -- for linear
 * acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian1DAcousticLinX(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);


/**
 * \brief Compute a 1D gaussian along y: A*exp(-(y-y0)^2/(2*var_y) + B -- for linear
 * acoustic.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = y_0,
 * coeffs[2] = var_y, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian1DAcousticLinY(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);


#endif // linAcoustic_boundaryCondition_hpp_included
