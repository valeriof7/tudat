/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120220    A. Ronse          File created.
 *      120307    K. Kumar          Updated description of circular restricted three-body problem.
 *      120321    K. Kumar          Moved circular restricted three-body problem inside gravitation
 *                                  namespace, which is inside a new astrodynamics namespace.
 *      120515    P. Musegaas       Added states namespace.
 *      120525    K. Kumar          Moved basic_mathematics namespace to tudat; added
 *                                  numerical_derivatives.
 *
 *    References
 */

#ifndef TUDAT_NAMESPACE_DESCRIPTIONS_H
#define TUDAT_NAMESPACE_DESCRIPTIONS_H

namespace tudat
{

//! Astrodynamics namespace.
/*!
 * The astrodynamics namespace is the top-level namespace for namespace related to astrodynamics
 * functionality.
 */
namespace astrodynamics
{

//! Acceleration models namespace.
/*!
 * In the acceleration models namespace, the various models for computing accelerations
 * (gravitational, aerodynamic, etc. ) are defined as free functions. In future, these will be
 * wrapped in classes too.
 */
namespace acceleration_models
{

} // namespace acceleration_models

//! Force models namespace.
/*!
 * In the force models namespace, the various models for computing forces
 * (gravitational, aerodynamic, etc. ) are defined as free functions. In future, these will be
 * wrapped in classes too.
 */
namespace force_models
{

} // namespace force_models

//! Gravitation namespace.
/*!
 * The gravitation namespace contains a base class for gravity field models. Derived classes to
 * generate a spherical harmonics (only interface implemented) and central gravity field are also
 * included. The namespace also includes functions to compute the gravitational force and
 * acceleration resulting from point masses. For the Circular Restricted Three-Body Problem
 * (CRTBP), this namespace also includes a function to compute Jacobi energy, a class to compute
 * the locations of the five Lagrange libration points, functions to convert from normalized,
 * dimensionless units to dimensional units, and a class to compute the state derivatives.
 */
namespace gravitation
{

//! Circular Restricted Three-Body Problem (CRTBP) namespace.
/*!
 * This namespace includes a function to compute Jacobi energy, a class to compute the locations of
 * the five Lagrange libration points, functions to convert from normalized, dimensionless state
 * and time to dimensional counterparts, and a class to compute the state derivatives.
 */
namespace circular_restricted_three_body_problem
{

} // namespace circular_restricted_three_body_problem

} // namespace gravitation

//! Moment models namespace.
/*!
 * In the moment models namespace, the various models for computing moments
 * (aerodynamic, due to general resultant force, etc. ) are defined as free functions. In future,
 * these will be wrapped in classes too.
 */
namespace moment_models
{

} // namespace moment_models

//! Rotational acceleration models namespace
/*!
 * In the rotational acceleration models namespace, the various models for computing
 * rotational accelerations (aerodynamic, etc. ) are defined as free functions. In future, these
 * will be wrapped in classes too.
 */
namespace rotational_acceleration_models
{

} // namespace rotational_acceleration_models

//! States namespace.
/*!
 * This namespace contains the state class and associated methods. Also this namespace contains
 * the Keplerian elements, Cartesian elements and Cartesian velocity and position elements.
 */
namespace states
{

} // namespace states

} // namespace astrodynamics

//! Aerodynamics namespace.
/*!
 * The aerodynamics namespace, containing constants and functions to compute flow characteristics
 * in super- and hypersonic flows, e.g., the Prandtl-Meyer function, various local surface
 * inclination methods, shock equations etc.
 */
namespace aerodynamics
{

} // namespace aerodynamics

//! Ephemerides namespace
/*!
 * The ephemerides namespace, containing the functionalities needed to obtain and propagate 
 * ephemerides. At this time, the JPL approximate planet positions are implemented as well as their
 * circular coplanar approximations.
 */
namespace ephemerides
{

} // namespace ephemerides

//! Orbital element conversions namespace.
/*!
 * The orbital element conversions namespace that contains free functions and classes to
 * convert between orbital elements. E.g., mean anomaly to eccentric anomaly, eccentric anomaly to
 * true anomaly etc.
 */
namespace orbital_element_conversions
{

} // namespace orbital_element_conversions

//! Reference frame namespace.
/*!
 * The reference frame namespace, containing functions to obtain the transformation quaternion or
   matrix between various reference frames. Includes rotating planetocentric, inertial,
   aerodynamic, body and local vertical reference frames.
 */
namespace reference_frame_transformations
{

} // namespace reference_frame_transformations

//! Output namespace.
/*!
 * Output namespace, containing functions to write state history to a file, get the root-directory
 * of the Tudat Package library, list all files in directory, read a file into a data matrix and
 * write surface geometry to a file. It
 * also contains a class containing TLE variables and a class to read TLE catalog files.
 */
namespace output
{

} // namespace output

//! Basic mathematics namespace.
/*!
 * The basic mathematics namespace, containing coordinate conversions, a nearest neighbor search
 * algorithm, and functions to compute numerical derivatives.
 */
namespace basic_mathematics
{

//! Coordinate conversions namespace.
/*!
 * This namespace contains functions to convert between different coordinates. At this point
 * conversion between cylindrical and Cartesian coordinates is included.
 */
namespace coordinate_conversions
{

} // namespace coordinate_conversions

//! numerical derivatives namespace.
/*!
 * This namespace contains functions to compute numerical derivatives. Currently, computations are
 * included to 2nd-, 4th- and 8th-order.
 */
namespace numerical_derivatives
{

} // namespace numerical_derivatives

} // namespace basic_mathematics

//! Mathematics namespace.
/*!
 * The mathematics namespace, containing various mathematical functions and related namespaces
 */
namespace mathematics
{

//! Interpolators namespace.
/*!
 * The interpolators namespace, containing the cubic spline interpolation class and functions to
 * perform a linear interpolation from data in map or vector form.
 * neighbour.
 */
namespace interpolators
{

} // namespace interpolators

//! Numerical integrators namespace.
/*!
 * The numerical integrators namespace, containing a class that implements the Euler integrator,
 * a struct that defines the coefficients of various Runge-Kutta integrators
 * ( RKF45, RKF56, RKF78, DOPRI87 ) and a class that implements the Runge-Kutta variable
 * step size integrator.
 */
namespace numerical_integrators
{

} // namespace interpolators

//! Statistics namespace.
/*!
 * The statistics namespace, containing functions to compute the average, standard deviation,
 * sample mean and sample variance.
 */
namespace statistics
{

} // namespace statistics

} // namespace mathematics

} // namespace tudat

#endif // TUDAT_NAMESPACE_DESCRIPTIONS_H