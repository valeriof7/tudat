/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H
#define TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating one-way range observables.
/*!
 *  Class for simulating one-way range, based on light-time and light-time corrections.
 *  The one-way range is defined as the light time multiplied by speed of light.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double,
          typename TimeType = double >
class OneWayRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    OneWayRangeObservationModel(
            const boost::shared_ptr< observation_models::LightTimeCalculator
            < ObservationScalarType, TimeType > > lightTimeCalculator,
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType >( oneWayRange, observationBiasCalculator ),
      lightTimeCalculator_( lightTimeCalculator ){ }

    //! Destructor
    ~OneWayRangeObservationModel( ){ }

    //! Function to compute ideal one-way range observation at given time.
    /*!
     *  This function compute ideal the one-way observation at a given time. The time argument can be either the reception
     *  or transmission time (defined by linkEndAssociatedWithTime input) Note that this observable does include e.g.
     *  light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Calculated observed one-way range value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        // Check link end associated with input time.
        bool isTimeAtReception = -1;
        if( linkEndAssociatedWithTime == receiver )
        {
            isTimeAtReception = 1;
        }
        else if( linkEndAssociatedWithTime == transmitter )
        {
            isTimeAtReception = 0;
        }
        else
        {
            throw std::runtime_error(
                        "Error when calculating one way range observation, link end is not transmitter or receiver" );
        }

        // Calculate light-time and multiply by speed of light in vacuum.
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) <<
                 lightTimeCalculator_->calculateLightTime( time, isTimeAtReception ) *
                 physical_constants::getSpeedOfLight< ObservationScalarType >( ) ).finished( );
    }

    //! Function to compute one-way range observable without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. Note that this observable does include light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        ObservationScalarType observation = TUDAT_NAN;
        TimeType transmissionTime = TUDAT_NAN, receptionTime = TUDAT_NAN;

        // Check link end associated with input time and compute observable
        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            observation = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, 1 );
            transmissionTime = time - observation;
            receptionTime = time;
            break;

        case transmitter:
            observation = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, 0 );
            transmissionTime = time;
            receptionTime = time + observation;
            break;
        default:
            std::string errorMessage = "Error, cannot have link end type: " +
                    boost::lexical_cast< std::string >( linkEndAssociatedWithTime ) + "for one-way range";
            throw std::runtime_error( errorMessage );
        }

        // Convert light time to range.
        observation *= physical_constants::getSpeedOfLight< ObservationScalarType >( );

        // Set link end states and times.
        linkEndTimes.push_back( static_cast< double >( transmissionTime ) );
        linkEndTimes.push_back( static_cast< double >( receptionTime ) );

        linkEndStates.push_back( transmitterState.template cast< double >( ) );
        linkEndStates.push_back( receiverState.template cast< double >( ) );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << observation ).finished( );
    }

    //! Function to get the object to calculate light time.
    /*!
     * Function to get the object to calculate light time.
     * \return Object to calculate light time.
     */
    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    lightTimeCalculator_;

    //! Pre-declared receiver state, to prevent many (de-)allocations
    StateType receiverState;

    //! Pre-declared transmitter state, to prevent many (de-)allocations
    StateType transmitterState;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H
