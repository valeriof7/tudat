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
 *      120203    B. Tong Minh      File created.
 *      120327    K. Kumar          Added Runge-Kutta 87 (Dormand and Prince) coefficients; added
 *                                  lower-, higher-order, and order to integrate variables.
 *
 *    References
 *      The Mathworks, Inc. RKF78, Symbolic Math Toolbox, 2012.
 *      Fehlberg, E. Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas With
 *          Stepsize Control, Marshall Spaceflight Center, NASA TR R-278, 1968.
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *
 *    Notes
 *      The naming of the coefficient sets follows (Montenbruck and Gill, 2005).
 *
 */

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace numerical_integrators
{

//! Initialize RKF45 coefficients.
void initializeRungeKuttaFehlberg45Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg45Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg45Coefficients.lowerOrder = 4;
    rungeKuttaFehlberg45Coefficients.higherOrder = 5;
    rungeKuttaFehlberg45Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 6, 5 );
    rungeKuttaFehlberg45Coefficients.aCoefficients( 1, 0 ) = 1.0 / 4.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 0 ) = 3.0 / 32.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 1 ) = 9.0 / 32.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 0 ) =  1932.0 / 2197.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 1 ) = -7200.0 / 2197.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 2 ) =  7296.0 / 2197.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 0 ) = 439.0 / 216.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 1 ) = -8.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 2 ) = 3680.0 / 513.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 3 ) = -845.0 / 4104.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 0 ) = -8.0 / 27.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 1 ) = 2.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 2 ) = -3544.0 / 2565.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 3 ) = 1859.0 / 4104.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 4 ) = -11.0 / 40.0;


    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.cCoefficients = Eigen::VectorXd::Zero( 6 );
    rungeKuttaFehlberg45Coefficients.cCoefficients( 0 ) = 0.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 1 ) = 1.0 / 4.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 2 ) = 3.0 / 8.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 3 ) = 12.0 / 13.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 4 ) = 1.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 5 ) = 1.0 / 2.0;


    // Define b-coefficients for the Runge-Kutta method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 6 );
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 0 ) = 25.0 / 216.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 2 ) = 1408.0 / 2565.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 3 ) = 2197.0 / 4104.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 4 ) = -1.0 / 5.0;

    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 0 ) = 16.0 / 135.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 2 ) = 6656.0 / 12825.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 3 ) = 28561.0 / 56430.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 4 ) = -9.0 / 50.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 5 ) = 2.0 / 55.0;
}

//! Initialize RKF56 coefficients.
void initializeRungeKuttaFehlberg56Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg56Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg56Coefficients.lowerOrder = 5;
    rungeKuttaFehlberg56Coefficients.higherOrder = 6;
    rungeKuttaFehlberg56Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 8, 7 );
    rungeKuttaFehlberg56Coefficients.aCoefficients( 1, 0 ) = 1.0 / 6.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 2, 0 ) =  4.0 / 75.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 2, 1 ) = 16.0 / 75.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 0 ) =  5.0 / 6.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 1 ) = -8.0 / 3.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 2 ) =  5.0 / 2.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 0 ) =  -8.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 1 ) = 144.0 / 25.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 2 ) =  -4.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 3 ) =  16.0 / 25.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 0 ) = 361.0 / 320.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 1 ) = -18.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 2 ) = 407.0 / 128.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 3 ) = -11.0 / 80.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 4 ) =  55.0 / 128.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 0 ) = -11.0 / 640.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 1 ) =   0.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 2 ) =  11.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 4 ) =  11.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 5 ) =   0.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 0 ) = 93.0 / 640.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 1 ) = -18.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 2 ) = 803.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 4 ) = 99.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 5 ) = 0.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 6 ) = 1.0;

    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.cCoefficients = Eigen::VectorXd::Zero( 8 );
    rungeKuttaFehlberg56Coefficients.cCoefficients( 1 ) = 1.0 / 6.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 2 ) = 4.0 / 15.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 3 ) = 2.0 / 3.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 4 ) = 4.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 5 ) = 1.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 7 ) = 1.0;

    // Define b-coefficients for the Runge-Kutta method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 8 );
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 0 ) = 31.0 / 384.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 2 ) = 1125.0 / 2816.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 3 ) = 9.0 / 32.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 4 ) = 125.0 / 768.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 5 ) = 5.0 / 66.0;

    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 0 ) = 7.0 / 1408.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 2 ) = 1125.0 / 2816.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 3 ) = 9.0 / 32.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 4 ) = 125.0 / 768.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 6 ) = 5.0 / 66.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 7 ) = 5.0 / 66.0;
}

//! Initialize RKF78 coefficients.
void initializeRungeKuttaFehlberg78Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg78Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg78Coefficients.lowerOrder = 7;
    rungeKuttaFehlberg78Coefficients.higherOrder = 8;
    rungeKuttaFehlberg78Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 7
    // with an embedded 8th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 13, 12 );
    rungeKuttaFehlberg78Coefficients.aCoefficients( 1, 0 ) = 2.0 / 27.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 0 ) = 1.0 / 36.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 1 ) = 1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 0 ) = 1.0 / 24.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 2 ) = 1.0 / 8.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 0 ) = 5.0 / 12.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 2 ) = -25.0 / 16.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 3 ) = 25.0 / 16.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 0 ) = 1.0 / 20.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 3 ) = 1.0 / 4.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 4 ) = 1.0 / 5.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 0 ) = -25.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 3 ) = 125.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 4 ) = -65.0 / 27.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 5 ) = 125.0 / 54.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 0 ) = 31.0 / 300.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 3 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 4 ) = 61.0 / 225.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 5 ) = -2.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 6 ) = 13.0 / 900.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 0 ) = 2.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 3 ) = -53.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 4 ) = 704.0 / 45.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 5 ) = -107.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 6 ) = 67.0 / 90.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 7 ) = 3.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 0 ) = -91.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 3 ) = 23.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 4 ) = -976.0 / 135.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 5 ) = 311.0 / 54.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 6 ) = -19.0 / 60.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 7 ) = 17.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 8 ) = -1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 0 ) = 2383.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 5 ) = -301.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 6 ) = 2133.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 7 ) = 45.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 8 ) = 45.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 9 ) = 18.0 / 41.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 0 ) = 3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 3 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 4 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 5 ) = -6.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 6 ) = -3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 7 ) = -3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 8 ) = 3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 9 ) = 6.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 10 ) = 0.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 0 ) = -1777.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 5 ) = -289.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 6 ) = 2193.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 7 ) = 51.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 8 ) = 33.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 9 ) = 12.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 10 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 11 ) = 1.0;

    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.cCoefficients = Eigen::VectorXd::Zero( 13 );
    rungeKuttaFehlberg78Coefficients.cCoefficients( 0 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 1 ) = 2.0 / 27.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 2 ) = 1.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 3 ) = 1.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 4 ) = 5.0 / 12.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 5 ) = 1.0 / 2.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 6 ) = 5.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 7 ) = 1.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 8 ) = 2.0 / 3.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 9 ) = 1.0 / 3.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 10 ) = 1.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 11 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 12 ) = 1.0;

    // Define b-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 13 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 0 ) = 41.0 / 840.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 5 ) = 34.0 / 105.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 6 ) = 9.0 / 35.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 7 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 6 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 8 ) = 9.0 / 280.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 9 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 8 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 10 ) = 41.0 / 840.0;

    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 5 ) = 34.0 / 105.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 6 ) = 9.0 / 35.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 7 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 6 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 8 ) = 9.0 / 280.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 9 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 8 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 11 ) = 41.0 / 840.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 12 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 11 );
}

//! Initialize RK87 (Dormand and Prince) coefficients.
void initializerungeKutta87DormandPrinceCoefficients(
        RungeKuttaCoefficients& rungeKutta87DormandPrinceCoefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta87DormandPrinceCoefficients.lowerOrder = 7;
    rungeKutta87DormandPrinceCoefficients.higherOrder = 8;
    rungeKutta87DormandPrinceCoefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::higher;

    // This coefficient set is taken from (Montenbruck and Gill, 2005).

    // a-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 13, 12 );

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 1, 0 ) = 1./18.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 2, 0 ) = 1./48.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 2, 1 ) = 1./16.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 3, 0 ) = 1./32.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 3, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 3, 2 ) = 3./32.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 0 ) = 5./16.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 2 ) = -75./64.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 3 ) = 75./64.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 0 ) = 3./80.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 3 ) = 3./16.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 4 ) = 3./20.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 0 ) = 29443841./614563906.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 3 ) = 77736538./692538347.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 4 ) = -28693883./1125000000.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 5 ) = 23124283./1800000000.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 0 ) = 16016141./946692911.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 3 ) = 61564180./158732637.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 4 ) = 22789713./633445777.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 5 ) = 545815736./2771057229.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 6 ) = -180193667./1043307555.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 0 ) = 39632708./573591083.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 3 ) = -433636366./683701615.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 4 ) = -421739975./2616292301.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 5 ) = 100302831./723423059.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 6 ) = 790204164./839813087.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 7 ) = 800635310./3783071287.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 0 ) = 246121993./1340847787.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 3 ) = -37695042795./15268766246.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 4 ) = -309121744./1061227803.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 5 ) = -12992083./490766935.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 6 ) = 6005943493./2108947869.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 7 ) = 393006217./1396673457.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 8 ) = 123872331./1001029789.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 0 ) = -1028468189./846180014.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 3 ) = 8478235783./508512852.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 4 ) = 1311729495./1432422823.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 5 ) = -10304129995./1701304382.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 6 ) = -48777925059./3047939560.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 7 ) = 15336726248./1032824649.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 8 ) = -45442868181./3398467696.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 9 ) = 3065993473./597172653.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 0 ) = 185892177./718116043.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 3 ) = -3185094517./667107341.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 4 ) = -477755414./1098053517.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 5 ) = -703635378./230739211.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 6 ) = 5731566787./1027545527.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 7 ) = 5232866602./850066563.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 8 ) = -4093664535./808688257.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 9 ) = 3962137247./1805957418.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 10 ) = 65686358./487910083.;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 0 ) = 403863854./491063109.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 3 ) = -5068492393./434740067.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 4 ) = -411421997./543043805.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 5 ) = 652783627./914296604.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 6 ) = 11173962825./925320556.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 7 ) = -13158990841./6184727034.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 8 ) = 3936647629./1978049680.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 9 ) = -160528059./685178525.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 10 ) = 248638103./1413531060.;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 11 ) = 0.;

    // c-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.cCoefficients = Eigen::VectorXd::Zero( 13 );

    rungeKutta87DormandPrinceCoefficients.cCoefficients( 0 ) = 0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 1 ) = 1./18.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 2 ) = 1./12.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 3 ) = 1./8.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 4 ) = 5./16.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 5 ) = 3./8.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 6 ) = 59./400.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 7 ) = 93./200.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 8 ) = 5490023248./9719169821. ;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 9 ) = 13./20.;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 10 ) = 1201146811./1299019798. ;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 11 ) = 1.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 12 ) = 1.0;

    // b-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 13 );

    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 0 ) = 14005451./335480064.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 3 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 4 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 5 ) = -59238493./1068277825.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 6 ) = 181606767./758867731.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 7 ) = 561292985./797845732.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 8 ) = -1041891430./1371343529.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 9 ) = 760417239./1151165299.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 10 ) = 118820643./751138087.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 11 ) = -528747749./2220607170.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 12 ) = 1./4.;

    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 0 ) = 13451932./455176623.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 1 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 2 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 3 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 4 ) = 0.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 5 ) = -808719846./976000145.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 6 ) = 1757004468./5645159321.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 7 ) = 656045339./265891186.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 8 ) = -3867574721./1518517206.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 9 ) = 465885868./322736535.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 10 ) = 53011238./667516719.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 11 ) = 2./45.;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 12 ) = 0.;
}

//! Get coefficients for a specified coefficient set
const RungeKuttaCoefficients& RungeKuttaCoefficients::get(
        RungeKuttaCoefficients::CoefficientSets coefficientSet )
{
    static RungeKuttaCoefficients rungeKuttaFehlberg45Coefficients,
            rungeKuttaFehlberg56Coefficients, rungeKuttaFehlberg78Coefficients,
            rungeKutta87DormandPrinceCoefficients;

    switch ( coefficientSet )
    {
    case rungeKuttaFehlberg45:
        if ( rungeKuttaFehlberg45Coefficients.higherOrder != 5 )
        {
            initializeRungeKuttaFehlberg45Coefficients( rungeKuttaFehlberg45Coefficients );
        }
        return rungeKuttaFehlberg45Coefficients;

    case rungeKuttaFehlberg56:
        if ( rungeKuttaFehlberg56Coefficients.higherOrder != 6 )
        {
            initializeRungeKuttaFehlberg56Coefficients( rungeKuttaFehlberg56Coefficients );
        }
        return rungeKuttaFehlberg56Coefficients;

    case rungeKuttaFehlberg78:
        if ( rungeKuttaFehlberg78Coefficients.higherOrder != 8 )
        {
            initializeRungeKuttaFehlberg78Coefficients( rungeKuttaFehlberg78Coefficients );
        }
        return rungeKuttaFehlberg78Coefficients;

    case rungeKutta87DormandPrince:
        if ( rungeKutta87DormandPrinceCoefficients.higherOrder != 8 )
        {
            initializerungeKutta87DormandPrinceCoefficients(
                        rungeKutta87DormandPrinceCoefficients );
        }
        return rungeKutta87DormandPrinceCoefficients;

    default: // The default case will never occur because CoefficientsSet is an enum.
        throw RungeKuttaCoefficients( );
    }
}

} // namespace numerical_integrators
} // namespace tudat
