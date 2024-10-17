#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>


#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"



namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_clock_parameter_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::statistics;
using namespace tudat::system_models;
using namespace tudat::propagators;
using namespace tudat::ground_stations;
using namespace tudat::basic_astrodynamics;

template< typename StateScalarType = double , typename TimeType = double >
Eigen::VectorXd  executeParameterEstimation( )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = TimeType( 1.2E7 );
    double maximumTimeStep = 3600.0;

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );

    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< AutoGeneratedTabulatedEphemerisSettings >(
            bodySettings.at( "Earth" )->ephemerisSettings, initialEphemerisTime - 4.0 * 86400.0, finalEphemerisTime + 4.0 * 86400, 3600.0 );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );

    createGroundStations( bodies, groundStationsToCreate );
    std::shared_ptr< GroundStation > grazStation = bodies.at( "Earth" )->getGroundStation( "Graz" );

    groundStationsToCreate.clear( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) << -4444.0,
                    unit_conversions::convertDegreesToRadians< double >( -4.49 ),
                    unit_conversions::convertDegreesToRadians< double >( 137.42 ) ).finished( );

    createGroundStations( bodies, groundStationsToCreate, coordinate_conversions::geodetic_position );
    std::shared_ptr< GroundStation > mslStation = bodies.at( "Mars" )->getGroundStation( "MSL" );

    // Create graz time error arcs.
    std::vector< Time > grazClockErrorArcTimes;
    double grazSingleArcLength = 0.5E6;
    Time currentTime = Time( initialEphemerisTime );
    while( currentTime < Time( finalEphemerisTime ) )
    {
        grazClockErrorArcTimes.push_back( currentTime );
        currentTime += grazSingleArcLength;
    }
    grazClockErrorArcTimes.push_back( finalEphemerisTime + grazSingleArcLength );
    grazClockErrorArcTimes.push_back( finalEphemerisTime + 2.0 * grazSingleArcLength );

    // Set errors for all arcs of ground station.
    //std::vector< double > grazArcPolynomialErrors;
    //grazArcPolynomialErrors.push_back( 1.0E-3 );
    //grazArcPolynomialErrors.push_back( -2.0E-9 );
    //grazArcPolynomialErrors.push_back( 3.1E-16 );

    std::vector< std::vector< double > > grazArcPolynomialErrors;
    for( unsigned int i = 0; i <5; i++ )
    {
        std::vector< double > currentErrors;
        currentErrors.push_back( ( static_cast< int >( i ) + 1.0 ) * 1.5 * 1.0E-3 );
        currentErrors.push_back( ( static_cast< int >( i ) + 1.0 ) * 1.5 * -2.0E-9 );
        currentErrors.push_back( ( static_cast< int >( i ) + 1.0 ) * 1.5 * 3.1E-16  );


        grazArcPolynomialErrors.push_back( currentErrors );
    }

    // Create and set timing system of graz.
    std::shared_ptr< TimingSystem > grazTimingSystem = std::make_shared< TimingSystem >( grazClockErrorArcTimes, grazArcPolynomialErrors );
    grazStation->setTimingSystem( grazTimingSystem );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToEstimate;
    bodiesToEstimate.push_back( "Earth" );
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralbodies;

    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";

        centralbodies[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >
                    ( 0.0, 10.0,
                      rungeKuttaFehlberg78, 1.0E-3, 1.0E3, 1.0E-12, 1.0E-12 );


    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                    propagators::getInitialStateOfBody( "Earth", "SSB", bodies, initialEphemerisTime ).template cast< StateScalarType >( ),
                    TimeType( initialEphemerisTime ),
                      integratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );

    // Set parameters that are to be estimated.
   std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< StateScalarType, TimeType >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< MultiArcPolynomialClockCorrectionsParameterSettings >
                              ( "Earth", "Graz",
                                std::vector< int >( { 1, 0, 2 } ),
                                std::vector< int >( { 0, 1, 2, 3 }) ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodies );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds2;
    linkEnds2.resize( 2 );
    linkEnds2[ 0 ][ transmitter ] = LinkEndId( "Earth", "Graz" );
    linkEnds2[ 0 ][ receiver ] = LinkEndId( "Mars", "MSL" );

    linkEnds2[ 1 ][ receiver ] = LinkEndId( "Earth", "Graz" );
    linkEnds2[ 1 ][ transmitter ] = LinkEndId( "Mars", "MSL" );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsMap;
    linkEndsMap[ one_way_range ] = linkEnds2;


    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
//    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "Mars", "MSL" ) );
    biasSettingsList.push_back( std::make_shared< TiminigSystemBiasSettings >( "Earth", "Graz" ) );

    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
            one_way_range, linkEnds2[ 0 ], nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList )  ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
            one_way_range, linkEnds2[ 1 ], nullptr, std::make_shared< MultipleObservationBiasSettings >( biasSettingsList )  ) );
    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                    bodies, parametersToEstimate,
                    observationSettingsList, propagatorSettings );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    double observationTimeStep = 1000.0;
    TimeType observationTime = Time( initialEphemerisTime + 10.0E4 );
    int numberOfObservations = 1800;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservations );

    for( int i = 0; i < numberOfObservations; i++ )
    {
        initialObservationTimes[ i ] = observationTime;
        observationTime += observationTimeStep;
    }

    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;

    initialObservationTimes = utilities::addScalarToVector( initialObservationTimes, 30.0 );
    measurementSimulationInput.push_back(
            std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                    one_way_range, linkEnds2[ 0 ], initialObservationTimes, receiver ) );


    singleObservableSimulationInput.clear( );

    std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > >  observationsAndTimes =
            simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
    }


    for( unsigned int i = bodiesToEstimate.size( ) * 6; i < initialParameterEstimate.rows( ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-1 );
    }


    parametersToEstimate->resetParameterValues( initialParameterEstimate );


    std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType> >(
                observationsAndTimes );

//    std::shared_ptr< orbit_determination::TimeConverter > timeConverter =
//            std::make_shared< orbit_determination::TimeConverter >(
//                grazStation, grazStation, local_proper_time_scale, local_proper_time_scale, bodies, 1, grazTimingSystem );
//
//    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< orbit_determination::TimeConverter > > > timeConverterMap;
//    timeConverterMap[ oneWayRange ][ linkEnds2[ 1 ] ] = timeConverter;

    std::shared_ptr< EstimationOutput< StateScalarType, TimeType > > podOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

    return ( podOutput->parameterEstimate_ - truthParameters ).cwiseQuotient( truthParameters ).template cast< double >( );
}

BOOST_AUTO_TEST_CASE( test_ClockParameterEstimation )
{
    Eigen::VectorXd relativeParameterError = executeParameterEstimation< long double, tudat::Time >( );
    for( unsigned int i = 6; i < relativeParameterError.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( relativeParameterError( i ) ), 1.0E-6 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

