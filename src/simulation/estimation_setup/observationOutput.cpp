/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/astro/observation_models/observationViabilityCalculator.h"

namespace tudat
{

namespace simulation_setup
{


void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings )
{
    if( isObservationDependentVariableGroundStationProperty( variableSettings ) )
    {
      
        std::string bodyName = variableSettings->relevantLinkEnd_.bodyName_;
        std::string stationName = variableSettings->relevantLinkEnd_.stationName_;

        if( bodies.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableId( variableSettings ) +
                                      ", could not find body " + bodyName );
        }

        if( bodies.at( bodyName )->getGroundStationMap( ).count( stationName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableId( variableSettings ) +
                                      ", could not find station " + stationName + " on body " + bodyName );
        }
    }
}


std::pair< int, int > getLinkEndStateTimeIndices(
    const observation_models::ObservableType observableType,
    const observation_models::LinkDefinition linkEnds,
    const observation_models::LinkEndId linkEndId,
    const observation_models::LinkEndType linkEndRole,
    const observation_models::LinkEndType originatingLinkEndRole,
    const IntegratedObservationPropertyHandling integratedObservableHandling )
{
    // Get link-end indices consistent with settings
    std::vector< std::pair< int, int > > currentStateTimeIndex =
        observation_models::getLinkStateAndTimeIndicesForLinkEnd(
            linkEnds.linkEnds_, observableType, linkEndId );

    // If no indices are found, throw exception
    if( currentStateTimeIndex.size( ) == 0 )
    {
        throw std::runtime_error( "Error in getting link and state ID for " +
                                  observation_models::getObservableName( observableType ) +
                                  " Could not find link end: (" + linkEndId.bodyName_ + ", " +
                                  linkEndId.stationName_ + ")" );
    }

    std::vector< std::pair< int, int > > stateTimeIndex;

    // Filter list of indices by link end role
    if( linkEndRole != observation_models::unidentified_link_end )
    {
        std::vector< int > linkEndIndices = observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
            observableType, linkEndRole, linkEnds.size( ) );
        for( unsigned int i = 0; i < currentStateTimeIndex.size( ); i++ )
        {
            for( unsigned int j = 0; j < linkEndIndices.size( ); j++ )
            {
                // If link end index and transmission index are equal, retain current link
                if( currentStateTimeIndex.at( i ).first == linkEndIndices.at( j ) )
                {
                    stateTimeIndex.push_back( currentStateTimeIndex.at( i ) );
                }
            }
        }
    }

    currentStateTimeIndex = stateTimeIndex;
    stateTimeIndex.clear( );

    // Filter list of indices by originating link end role
    if( originatingLinkEndRole != observation_models::unidentified_link_end )
    {
        std::vector< int > linkEndIndices = observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
            observableType, originatingLinkEndRole, linkEnds.size( ) );
        for( unsigned int i = 0; i < currentStateTimeIndex.size( ); i++ )
        {
            for( unsigned int j = 0; j < linkEndIndices.size( ); j++ )
            {
                // If originating link end index and transmission index are equal, retain current link
                if( currentStateTimeIndex.at( i ).second == linkEndIndices.at( j ) )
                {
                    stateTimeIndex.push_back( currentStateTimeIndex.at( i ) );
                }
            }
        }
    }

    currentStateTimeIndex = stateTimeIndex;
    stateTimeIndex.clear( );

    // Check if observation is integrated
    if( observation_models::isObservableOfIntegratedType( observableType ) )
    {
        if( currentStateTimeIndex.size( ) != 2 )
        {
            throw std::runtime_error( "Error when getting integrated observable state and time indices; 2 remaining indices required" );
        }
        else if( integratedObservableHandling == interval_undefined )
        {
            throw std::runtime_error( "Error when getting integrated observable state and time indices; choice of start or end of arc undefined" );
        }
        else if( integratedObservableHandling == interval_start )
        {
            stateTimeIndex.push_back( currentStateTimeIndex.at( 0 ) );
        }
        else if( integratedObservableHandling == interval_end )
        {
            stateTimeIndex.push_back( currentStateTimeIndex.at( 1 ) );
        }
    }
    else
    {
        if( currentStateTimeIndex.size( ) != 1 )
        {
            throw std::runtime_error( "Error when getting observable state and time indices; 1 index required" );
        }

        stateTimeIndex = currentStateTimeIndex;
    }

    return stateTimeIndex.at( 0 );
}

ObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds )
{
    // Check if relevant properties of environment exist
    checkObservationDependentVariableEnvironment( bodies, variableSettings );

    // Retrieve link-end ID for station
    std::string bodyName = variableSettings->relevantLinkEnd_.bodyName_;
    std::string stationName = variableSettings->relevantLinkEnd_.stationName_;

    std::pair< int, int > linkEndIndicesToUse = getLinkEndStateTimeIndices(
        observableType, linkEnds, variableSettings->relevantLinkEnd_,
        variableSettings->linkEndRole_, variableSettings->originatingLinkEndRole_,
        variableSettings->integratedObservableHandling_ );

    // Retrieve pointing angles calculator for station, and setup output function
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            bodies.at( bodyName )->getGroundStationMap( ).at( stationName )->getPointingAnglesCalculator( );
    ObservationDependentVariableFunction outputFunction;

    if( variableSettings->variableType_ == station_elevation_angle )
    {
        outputFunction = [=]( const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observationValue,
                const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySimulationSettings )
        {
            return ( Eigen::VectorXd( 1 ) << ground_stations::calculateGroundStationElevationAngle(
                        pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ) ).finished( );
        };
    }
    else if( variableSettings->variableType_ == station_azimuth_angle )
    {
        outputFunction = [=]( const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observationValue,
                const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySimulationSettings )
        {
            return ( Eigen::VectorXd( 1 ) << ground_stations::calculateGroundStationAzimuthAngle(
                        pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ) ).finished( );
        };
    }
    return outputFunction;
}


ObservationDependentVariableFunction getInterlinkObservationVariableFunction(
    const SystemOfBodies& bodies,
    const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings,
    const observation_models::ObservableType observableType,
    const observation_models::LinkDefinition linkEnds )
{
    std::pair< int, int > linkEndIndicesToUse = getLinkEndStateTimeIndices(
        observableType, linkEnds, linkEnds.at( variableSettings->startLinkEnd_ ),
        variableSettings->startLinkEnd_, variableSettings->endLinkEnd_,
        variableSettings->integratedObservableHandling_ );

    ObservationDependentVariableFunction outputFunction;

    switch( variableSettings->variableType_ )
    {
    case target_range:
    {
        if ( variableSettings->relativeBody_ != "" )
        {
            throw std::runtime_error(
                "Error when parsing target range observation dependent variable, relative body must not be defined" );
        }
        outputFunction = [ = ]( const std::vector<double> &linkEndTimes,
                                const std::vector<Eigen::Matrix<double, 6, 1> > &linkEndStates,
                                const Eigen::VectorXd &observationValue,
                                const std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> ancilliarySimulationSettings )
        {
            return ( Eigen::VectorXd( 1 ) << ( linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ) -
                                               linkEndStates.at( linkEndIndicesToUse.second ).segment( 0,
                                                                                                       3 )).norm( )).finished( );
        };
        break;
    }
    case body_avoidance_angle_variable:
    {
        if ( bodies.count( variableSettings->relativeBody_ ) != 0 )
        {
            throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                      variableSettings->relativeBody_ + ", body is not defined" );
        }

        if ( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
        {
            throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                      variableSettings->relativeBody_ + ", body has no ephemeris" );
        }

        outputFunction = [ = ]( const std::vector<double> &linkEndTimes,
                                const std::vector<Eigen::Matrix<double, 6, 1> > &linkEndStates,
                                const Eigen::VectorXd &observationValue,
                                const std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> ancilliarySimulationSettings )
        {
            double cosineOfAvoidanceAngle = observation_models::computeCosineBodyAvoidanceAngle(
                linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ),
                linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ),
                bodies.at( variableSettings->relativeBody_ )->getStateInBaseFrameFromEphemeris<double, double>(
                    ( linkEndIndicesToUse.first + linkEndIndicesToUse.second ) / 2.0 ).segment( 0, 3 ));
            if ( std::fabs( cosineOfAvoidanceAngle ) >= 1.0 )
            {
                return ( Eigen::Vector1d( )
                    << utilities::sgn( cosineOfAvoidanceAngle ) * mathematical_constants::PI / 2.0 ).finished( );
            }
            else
            {
                return ( Eigen::Vector1d( ) << std::acos( cosineOfAvoidanceAngle )).finished( );
            }
        };
        break;
    }
    case link_body_center_distance:
    {
        if ( bodies.count( variableSettings->relativeBody_ ) != 0 )
        {
            throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                      variableSettings->relativeBody_ + ", body is not defined" );
        }

        if ( bodies.at( variableSettings->relativeBody_ )->getEphemeris( ) == nullptr )
        {
            throw std::runtime_error( "Error when parsing body avoidance observation dependent variable w.r.t. " +
                                      variableSettings->relativeBody_ + ", body has no ephemeris" );
        }

        outputFunction = [ = ]( const std::vector<double> &linkEndTimes,
                                const std::vector<Eigen::Matrix<double, 6, 1> > &linkEndStates,
                                const Eigen::VectorXd &observationValue,
                                const std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> ancilliarySimulationSettings )
        {
            double minimumDistance = observation_models::computeMinimumLinkDistanceToPoint(
                linkEndStates.at( linkEndIndicesToUse.first ).segment( 0, 3 ),
                linkEndStates.at( linkEndIndicesToUse.second ).segment( 0, 3 ),
                bodies.at( variableSettings->relativeBody_ )->getStateInBaseFrameFromEphemeris<double, double>(
                    ( linkEndIndicesToUse.first + linkEndIndicesToUse.second ) / 2.0 ).segment( 0, 3 ));
            return ( Eigen::Vector1d( ) << minimumDistance ).finished( );
        };
        break;
    }
    case link_limb_distance:
    {
        std::shared_ptr< InterlinkObservationDependentVariableSettings > bodyCenterDistanceSettings =
            std::make_shared< InterlinkObservationDependentVariableSettings >(
                link_body_center_distance,
                variableSettings->startLinkEnd_, variableSettings->endLinkEnd_,
                variableSettings->integratedObservableHandling_, variableSettings->relativeBody_ );
        ObservationDependentVariableFunction linkCenterFunction = getInterlinkObservationVariableFunction(
            bodies, variableSettings, observableType, linkEnds );

        if ( bodies.at( variableSettings->relativeBody_ )->getShapeModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when parsing body link limb distance observation dependent variable w.r.t. " +
                                      variableSettings->relativeBody_ + ", body has no shape model" );
        }
        auto shapeModel = bodies.at( variableSettings->relativeBody_ )->getShapeModel( );
        outputFunction = [ = ]( const std::vector<double> &linkEndTimes,
                                const std::vector<Eigen::Matrix<double, 6, 1> > &linkEndStates,
                                const Eigen::VectorXd &observationValue,
                                const std::shared_ptr<observation_models::ObservationAncilliarySimulationSettings> ancilliarySimulationSettings )
        {
            return linkCenterFunction( linkEndTimes, linkEndStates, observationValue, ancilliarySimulationSettings ) -
                ( Eigen::Vector1d( ) << shapeModel->getAverageRadius( ) ).finished( );
        };
        break;
    }
    default:
        throw std::runtime_error( "Error when parsing interlink observation dependent variable, did not recognize variable" +
                              getObservationDependentVariableId( variableSettings ) );
    }
    return outputFunction;
}

ObservationDependentVariableFunction getObservationVectorDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkDefinition linkEnds )
{
    ObservationDependentVariableFunction outputFunction;
    switch( variableSettings->variableType_ )
    {
    case station_elevation_angle:
    {
        std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
            std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                variableSettings );
        if( angleSettings == nullptr )
        {
            throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_elevation_angle" );
        }
        outputFunction = getStationObservationAngleFunction(
            bodies, angleSettings, observableType, linkEnds );
        break;
    }
    case station_azimuth_angle:
    {
        std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
            std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                variableSettings );
        if( angleSettings == nullptr )
        {
            throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_azimuth_angle" );
        }

        outputFunction = getStationObservationAngleFunction(
            bodies, angleSettings, observableType, linkEnds );
        break;
    }
    case target_range:
    {
        std::shared_ptr< InterlinkObservationDependentVariableSettings > linkSettings =
            std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings );
        if( linkSettings == nullptr )
        {
            throw std::runtime_error( "Error in observation dependent variables, incorrect type found for target_range" );
        }

        outputFunction = getInterlinkObservationVariableFunction(
            bodies, linkSettings, observableType, linkEnds );
        break;
    }
    default:
        throw std::runtime_error( "Error when parsing vector observation dependent variable, did not recognize variable" +
                                  getObservationDependentVariableId( variableSettings ) );
    }
    return outputFunction;
}

void ObservationDependentVariableCalculator::addDependentVariable(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const SystemOfBodies& bodies )
{
    // Check if the requested dependent variable can be used for given link
//    if( doesObservationDependentVariableExistForGivenLink(
//                observableType_, linkEnds_.linkEnds_, variableSettings ) )
    {
        // Retrieve the current index in list of dependent variables and size of new parameter
        int currentIndex = totalDependentVariableSize_;
        int parameterSize = getObservationDependentVariableSize( variableSettings );

        // Create function to compute dependent variable
        ObservationDependentVariableFunction observationDependentVariableFunction =
                getObservationVectorDependentVariableFunction(
                    bodies, variableSettings, observableType_, linkEnds_ );

        // Create function to compute dependent variable and add to existing list
        ObservationDependentVariableAddFunction dependentVariableAddFunction = [=](
                Eigen::VectorXd& dependentVariables,
                const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observable,
                const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySimulationSettings )
        {
            // Check if computation is not overriding existing values
            for( int i = 0; i < parameterSize; i++ )
            {
                if( dependentVariables( currentIndex + i ) == dependentVariables( currentIndex + i ) )
                {
                    throw std::runtime_error( "Error when saving observation dependent variables; overriding existing value" );
                }
            }
            dependentVariables.segment( currentIndex, parameterSize ) = observationDependentVariableFunction(
                        linkEndTimes, linkEndStates, observable, ancilliarySimulationSettings );
        };
        
        // Add new dependent variable function and settings to list
        dependentVariableAddFunctions_.push_back( dependentVariableAddFunction );
        dependentVariableStartIndices_.push_back( totalDependentVariableSize_ );
        dependentVariableSizes_.push_back( parameterSize );
        settingsList_.push_back( variableSettings );
        totalDependentVariableSize_ += parameterSize;
    }
}

void ObservationDependentVariableCalculator::addDependentVariables(
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsList,
        const SystemOfBodies& bodies )
{
    for( unsigned int i = 0; i < settingsList.size( ); i++ )
    {
        addDependentVariable( settingsList.at( i ), bodies );
    }
}

std::pair< int, int > ObservationDependentVariableCalculator::getDependentVariableIndices(
        const std::shared_ptr< ObservationDependentVariableSettings > dependentVariables )
{
    std::pair< int, int > startAndSizePair = std::make_pair( 0, 0 );
    for( unsigned int i = 0; i < settingsList_.size( ); i++ )
    {
        if( getObservationDependentVariableId( settingsList_.at( i ) ) ==
                getObservationDependentVariableId( dependentVariables ) )
        {
            if( startAndSizePair.second == 0 )
            {
                startAndSizePair = std::make_pair(
                            dependentVariableStartIndices_.at( i ), dependentVariableSizes_.at( i ) );
            }
            else
            {
                throw std::runtime_error( "Error when finding observation dependent variable; multiple candidates found" );
            }
        }
    }
    return startAndSizePair;
}



Eigen::VectorXd ObservationDependentVariableCalculator::calculateDependentVariables(
        const std::vector< double >& linkEndTimes,
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const Eigen::VectorXd& observation,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > observationAncilliarySimulationSettings )
{
    Eigen::VectorXd dependentVariables = Eigen::VectorXd::Constant( totalDependentVariableSize_, TUDAT_NAN );

    for( unsigned int i = 0; i < dependentVariableAddFunctions_.size( ); i++ )
    {
        dependentVariableAddFunctions_.at( i )(
                    dependentVariables, linkEndTimes, linkEndStates, observation, observationAncilliarySimulationSettings );

    }
    return dependentVariables;
}




}


}
