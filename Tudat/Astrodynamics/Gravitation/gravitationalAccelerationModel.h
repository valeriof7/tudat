/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120209    K. Kumar          File created.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 */

#ifndef TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
#define TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H

#include <Eigen/Core>

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Compute gravitational acceleration.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with body2.
 * The basis for this gravitational acceleration is that body2 is a point mass,
 * generating acceleration due to Newton's gravitational force:
 * \f[
 *      \bar{a}_{gravity} = -\frac{G * m_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$G\f$ is the universal gravitational constant, \f$m_{2}\f$ is the mass of body2,
 * and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to body 2, with respect to an
 * inertial (barycentric) reference frame.
 * \param universalGravitationalParameter Universal gravitational constant [m^3 kg^-1 s^-2].
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to acceleration
 *          (body1) [m].
 * \param massOfBodyExertingAcceleration Mass of body exerting acceleration (body2) [kg].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 [m s^-2].
 * \sa computeGravitationalForce.
 */
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalParameter,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

//! Compute gravitational acceleration.
/*!
 * Computes gravitational acceleration experienced by body1, due to its interaction with another
 * body (body2). The basis for this gravitational acceleration is that the body exerting the
 * acceleration is a point mass, generating acceleration due to Newton's gravitational force:
 * \f[
 *      \bar{a}_{gravity} = -\frac{\mu_{2}}{r_{12}^{3}} * \bar{r}_{12}
 * \f]
 * where \f$\mu_{2}\f$ is the gravitational parameter of the body exerting acceleration,
 * and \f$\bar{r}_{12}\f$ is the relative position vector from body1 to body 2. with respect to an
 * inertial (barycentric) reference frame.
 * \param positionOfBodySubjectToAcceleration Position vector of body subject to acceleration
 *          (body1) [m].
 * \param gravitationalParameterOfBodyExertingAcceleration Gravitational parameter of body exerting
 *          acceleration (body2) [m^3 s^-2].
 * \param positionOfBodyExertingAcceleration Position vector of body exerting acceleration
 *          (body2) [m].
 * \return Gravitational acceleration exerted on body1 [m s^-2].
 * \sa computeGravitationalForce.
 */
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration );

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITATIONAL_ACCELERATION_MODEL_H
