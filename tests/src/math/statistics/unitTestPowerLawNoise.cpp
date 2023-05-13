
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include <ctime>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>

#include "tudat/math/statistics/powerLawNoiseGeneration.h"
#include "tudat/math/statistics/allanVariance.h"
#include "tudat/math/interpolators/linearInterpolator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_power_law_noise_generator )

using namespace tudat::statistics;
//BOOST_AUTO_TEST_CASE( testPowerLawNoise )
//{
//    time_t tstart, tend;
//    tstart = time(0);

//    std::vector< double > timeDomainNoise = generateScaleIndependentPowerLawNoise( -1.99, std::pow( 2, 12 ) + 1398373 );
//    output::writeVectorToFile( timeDomainNoise, "timeDomainNoise.dat"  );

//    tend = time(0);
//    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;
//}

//BOOST_AUTO_TEST_CASE( testInterpolatedPowerLawNoise )
//{
//    const double startTime = 0.0;
//    const double endTime = 1.0E7;
//    const double timeStep = 1.0;
//    const double power = 1.0;
//    const double amplitudeAtUnity = 1.0;

//    std::map< double, double > powerLawNoise = generatePowerLawNoise(
//                startTime, endTime, timeStep, power, amplitudeAtUnity );

//    std::vector< std::complex< double > > frequencyDomainPowerLawNoise =
//            fftw_interface::performFftOfRealData( utilities::createVectorFromMapValues( powerLawNoise ) );
//    //output::writeComplexVectorToFile( frequencyDomainPowerLawNoise,
//    //                                  "generatedNoise.dat" );

//    std::cout<<"T1"<<std::endl;

//    using namespace tudat::interpolators;

//    LinearInterpolatorDouble noiseInterpolator = LinearInterpolatorDouble( powerLawNoise );
//    double sampleTime = 1.0;

//    std::vector< double > overSampledNoise;
//    double currentTime = startTime;

//    while( currentTime < endTime - 1 )
//    {
//        overSampledNoise.push_back( noiseInterpolator.interpolate( currentTime + 0.2 ) );
//        currentTime += sampleTime;
//    }

//    std::vector< std::complex< double > > sampledNoiseInFrequencyDomain = fftw_interface::performFftOfRealData(
//                overSampledNoise );

//    //output::writeComplexVectorToFile( sampledNoiseInFrequencyDomain, "overSampledNoise.dat" );
//}

//BOOST_AUTO_TEST_CASE( testInterpolatedPowerLawNoise )
//{

//    std::vector< double > timeDomainNoise1 = generatePowerLawNoise( 1.0, std::pow( 2, 16 ), -1.0, 1.0 );
//    output::writeVectorToFile( timeDomainNoise1, "timeDomainNoise1.dat" );
//    std::vector< std::complex< double > > frequencyDomainNoise1 = fftw_interface::performFftOfRealData( timeDomainNoise1 );
//    output::writeComplexVectorToFile( frequencyDomainNoise1, "frequencyDomainNoise1.dat" );

//    std::vector< double > timeDomainNoise2 = generatePowerLawNoise( 1.0, std::pow( 2, 16 ), 2.0, 1.0 );
//    output::writeVectorToFile( timeDomainNoise2, "timeDomainNoise2.dat" );
//    std::vector< std::complex< double > > frequencyDomainNoise2 = fftw_interface::performFftOfRealData( timeDomainNoise2 );
//    output::writeComplexVectorToFile( frequencyDomainNoise2, "frequencyDomainNoise2.dat" );


//}

double getExpectedAmplitudeAtUnitTime( const int power, const double maximumFrequency )
{
    double expectedAmplitde = 0.0;
    switch( power )
    {
    case -4:
        expectedAmplitde = 2.0 * mathematical_constants::PI * mathematical_constants::PI / 3.0;
        break;
    case -3:
        expectedAmplitde = 2.0 * std::log( 2.0 );
        break;
    case -2:
        expectedAmplitde = 0.5;
        break;
    case 0:
        expectedAmplitde = 3.0 * maximumFrequency / ( 4.0 * mathematical_constants::PI * mathematical_constants::PI );
        break;
    default:
        std::cerr<<"Error when getting expected allan variance amplitude at unit time interval, found power "<<power<<std::endl;
        break;
    }
    return expectedAmplitde;
}

double getExpectedAllanVarianceExponent( const int phasePower )
{
    int allanVarianceExponent = 0.0;
    switch( phasePower )
    {
    case -4:
        allanVarianceExponent = 1;
        break;
    case -3:
        allanVarianceExponent = 0;
        break;
    case -2:
        allanVarianceExponent = -1;
        break;
    case -1:
        allanVarianceExponent = -2;
        break;
    case 0:
        allanVarianceExponent = -2;
        break;
    default:
        std::cerr<<"Error when getting expected allan variance exponent, found power "<<phasePower<<std::endl;
        break;
    }
    return allanVarianceExponent;
}

BOOST_AUTO_TEST_CASE( testPowerLawNoiseSimulation )
{
    std::vector< double > maximumFrequencies = boost::assign::list_of( 128.0 );
    int counter = 0;
    double timeIntervalRatio = 0.0;
    double allanVarianceRatio = 0.0;
    std::map< double, double >::iterator startIterator;
    std::map< double, double >::iterator endIterator;

    for( int i = -4; i <=0; i++ )
    {
        if( i!= -1 )
        {
            for( unsigned int j = 0; j < maximumFrequencies.size( ); j++ )
            {
                unsigned int numberOfRuns = 200;
                double currentValue = 0.0;
                double summedValue = 0.0;

                double currentExponentValue = 0.0;
                double summedExponentValue = 0.0;

                for( unsigned k = 0; k < numberOfRuns; k++ )
                {
                    std::pair< std::vector< double >, double > timeDomainNoiseOutput = generatePowerLawNoise(
                                maximumFrequencies.at( j ), maximumFrequencies.at( j ) / static_cast< double >( std::pow( 2.0, 16.0 ) ),
                                i, 1.0 / std::pow( 2.0 * mathematical_constants::PI , 2.0 ), static_cast< double >( counter ) );
                    std::map< double, double > allanVariance = statistics::calculateAllanVarianceOfTimeDataSet(
                                timeDomainNoiseOutput.first, timeDomainNoiseOutput.second );

                    startIterator = allanVariance.begin( );
                    endIterator = allanVariance.end( );
                    for( int m = 0; m < 4; m++ )
                    {
                        endIterator--;
                    }

                    timeIntervalRatio = endIterator->first / startIterator->first;
                    allanVarianceRatio = endIterator->second / startIterator->second;

                    currentExponentValue = std::log( allanVarianceRatio ) / std::log( timeIntervalRatio );
                    summedExponentValue += currentExponentValue;

                    interpolators::LinearInterpolator< double, double > allanVarianceInterpolator =
                            interpolators::LinearInterpolator< double, double >( allanVariance );
                    currentValue = allanVarianceInterpolator.interpolate( 1.0 );
                    //std::cout<<currentValue<<" ";

                    summedValue += currentValue;
                    //output::writeVectorToFile( timeDomainNoiseOutput.first, "timeDomainNoise" + boost::lexical_cast< std::string >( i ) + "_" +
                    //                              boost::lexical_cast< std::string >( k )+ ".dat" );
                    //output::writeDoubleMapToFile( allanVariance, "allanVariance_" + boost::lexical_cast< std::string >( i ) + "_" +
                    //                              boost::lexical_cast< std::string >( k )+ ".dat" );
                    counter++;
                    //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                }

                double expectedAmplitde = getExpectedAmplitudeAtUnitTime( i, maximumFrequencies.at( j ) );
                std::cout<<i<<" "<<j<<std::endl;
                BOOST_CHECK_CLOSE_FRACTION( expectedAmplitde, summedValue / static_cast< double >( numberOfRuns ), 0.02 );
                BOOST_CHECK_SMALL( std::fabs( summedExponentValue / static_cast< double >( numberOfRuns )  -
                                              static_cast< double >( getExpectedAllanVarianceExponent( i ) ) ), 0.06 );

            }
        }
    }
}

//BOOST_AUTO_TEST_CASE( testPowerLawNoiseSimulation )
//{
//    std::vector< double > maximumFrequencies = boost::assign::list_of( 4.0 );
//    int counter = 0;
//    double timeIntervalRatio = 0.0;
//    double allanVarianceRatio = 0.0;
//    std::map< double, double >::iterator startIterator;
//    std::map< double, double >::iterator endIterator;

//    for( int i = -4; i <=0; i++ )
//    {
//        if( i!= -1 )
//        {
//            for( unsigned int j = 0; j < maximumFrequencies.size( ); j++ )
//            {
//                unsigned int numberOfRuns = 1;
//                double currentValue = 0.0;
//                double summedValue = 0.0;

//                double currentExponentValue = 0.0;
//                double summedExponentValue = 0.0;

//                for( unsigned k = 0; k < numberOfRuns; k++ )
//                {
//                    std::pair< std::vector< double >, double > timeDomainNoiseOutput = generatePowerLawNoise(
//                                maximumFrequencies.at( j ), maximumFrequencies.at( j ) / static_cast< double >( std::pow( 2.0, 18.0 ) ),
//                                i, 1.0 / std::pow( 2.0 * mathematical_constants::PI , 2.0 ) );
//                    std::map< double, double > allanVariance = statistics::calculateAllanVarianceOfTimeDataSet(
//                                timeDomainNoiseOutput.first, timeDomainNoiseOutput.second );

//                    startIterator = allanVariance.begin( );
//                    endIterator = allanVariance.end( );
//                    for( int m = 0; m < 4; m++ )
//                    {
//                        endIterator--;
//                    }

//                    timeIntervalRatio = endIterator->first / startIterator->first;
//                    allanVarianceRatio = endIterator->second / startIterator->second;

//                    currentExponentValue = std::log( allanVarianceRatio ) / std::log( timeIntervalRatio );
//                    summedExponentValue += currentExponentValue;

//                    interpolators::LinearInterpolator< double, double > allanVarianceInterpolator =
//                            interpolators::LinearInterpolator< double, double >( allanVariance );
//                    currentValue = allanVarianceInterpolator.interpolate( 1.0 );
//                    //std::cout<<currentValue<<" ";
//                    std::cout<<timeDomainNoiseOutput.second<<std::endl;
//                    summedValue += currentValue;
//                    output::writeVectorToFile( timeDomainNoiseOutput.first, "timeDomainNoise" + boost::lexical_cast< std::string >( i ) + "_" +
//                                                  boost::lexical_cast< std::string >( k )+ ".dat" );
//                    //output::writeDoubleMapToFile( allanVariance, "allanVariance_" + boost::lexical_cast< std::string >( i ) + "_" +
//                    //                              boost::lexical_cast< std::string >( k )+ ".dat" );
//                    counter++;
//                    //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
//                }

//                double expectedAmplitde = getExpectedAmplitudeAtUnitTime( i, maximumFrequencies.at( j ) );
//                std::cout<<i<<" "<<j<<std::endl;
//                BOOST_CHECK_CLOSE_FRACTION( expectedAmplitde, summedValue / static_cast< double >( numberOfRuns ), 0.02 );
//                BOOST_CHECK_SMALL( std::fabs( summedExponentValue / static_cast< double >( numberOfRuns )  -
//                                              static_cast< double >( getExpectedAllanVarianceExponent( i ) ) ), 0.06 );

//            }
//        }
//    }
//}


BOOST_AUTO_TEST_SUITE_END( )

}

}
