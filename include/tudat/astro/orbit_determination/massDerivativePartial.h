/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MASSRATEPARTIAL_H
#define TUDAT_MASSRATEPARTIAL_H

#include "tudat/astro/orbit_determination/stateDerivativePartial.h"


namespace tudat
{

namespace orbit_determination
{

class MassRatePartial: public StateDerivativePartial
{
public:
    MassRatePartial( const std::string& body,
                     const basic_astrodynamics::AvailableMassRateModels massRateType ):
        StateDerivativePartial( propagators::body_mass_state, std::make_pair( body, "" ) ),
        body_( body ), massRateType_( massRateType )
    { }

    //! Virtual destructor.
    virtual ~MassRatePartial( ) { }

    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
    getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        // Initialize to empty function; 0 parameter size.
        std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
                partialFunction = std::make_pair( std::function< void( Eigen::Block< Eigen::MatrixXd > ) >( ), 0 );

        // Check if state dependency exists
        switch( integratedStateType )
        {
        case propagators::translational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting mass rate derivative model, cannot have reference point on body for translational dynamics" );
            }
            else if( body_ == stateReferencePoint .first && isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                         stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &MassRatePartial::wrtTranslationalStateOfBody, this, std::placeholders::_1 ), 6 );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                         stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &MassRatePartial::wrtTranslationalStateOfAdditionalBody, this, std::placeholders::_1, stateReferencePoint.first ), 6 );
            }
            break;
        }
        case propagators::rotational_state:
        {
            throw std::runtime_error( "Error when getting mass rate partial, rotational dynamics partial no yet implemented" );
            break;
        }
        case propagators::body_mass_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for body mass" );
            }
            else if( isMassRatePartialWrtMassNonZero( ) )
            {
                partialFunction = std::make_pair( std::bind( &MassRatePartial::wrtMassOfBody,
                                                             this, std::placeholders::_1 ), 1 );
            }
            break;
        }
        case propagators::custom_state:
        {
            break;
        }
        default:
            std::string errorMessage =
                    "Error when getting state derivative partial acceleration model, dynamics type " +
                    std::to_string( integratedStateType ) + "not recognized" ;
            throw std::runtime_error( errorMessage );
            break;
        }


        return partialFunction;
    }

    virtual void wrtMassOfBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {

    }

    virtual void wrtTranslationalStateOfBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {

    }

    virtual void wrtTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::string& bodyName )
    {

    }

    virtual bool isMassRatePartialWrtMassNonZero( ) = 0;

    std::string getBody( ) { return body_; }


    basic_astrodynamics::AvailableMassRateModels getMassRateType( )
    {
        return massRateType_;
    }

protected:

    std::string body_;

    basic_astrodynamics::AvailableMassRateModels massRateType_;
};

class FromThrustMassRatePartial: public MassRatePartial
{
public:
    FromThrustMassRatePartial( const std::string& body,
                               const std::shared_ptr< propulsion::FromThrustMassRateModel > massRateModel ):
        MassRatePartial( body, basic_astrodynamics::from_thrust_mass_rate_model )
    {
        thrustAccelerations_ = massRateModel->getThrustAccelerations( );

        for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
        {
            std::vector< std::shared_ptr< system_models::EngineModel > > thrustSources =
                    thrustAccelerations_.at( i )->getThrustSources( );
            std::vector< std::shared_ptr< system_models::EngineModel > > currentMassDependentThrustSources;
            for( unsigned int j = 0; j < thrustSources.size( ); j++ )
            {
                if( !thrustSources.at( i )->getThrustMagnitudeWrapper( )->modelIsForceBased( ) )
                {
                    currentMassDependentThrustSources.push_back( thrustSources.at( i ) );
                }
            }
            if( currentMassDependentThrustSources.size( ) > 0 )
            {
                accelerationBasedThrustSources_[ i ] = currentMassDependentThrustSources;

            }
        }
    }

    bool isMassRatePartialWrtMassNonZero( )
    {
        if( accelerationBasedThrustSources_.size( ) == 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    void wrtMassOfBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        for( auto it : accelerationBasedThrustSources_ )
        {
            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                partialMatrix( 0, 0 ) += it.second.at( i )->getCurrentMassRate( ) /
                        thrustAccelerations_.at( i )->getCurrentBodyMass( );
            }
        }
    }

    void wrtEngineSpecificImpulse(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::string& engineName )
    {
        std::shared_ptr< system_models::EngineModel > engineModel = engineModelList_.at( engineName );
        partialMatrix( 0, 0 ) += -engineModel->getCurrentMassRate( ) /
                engineModel->getThrustMagnitudeWrapper( )->getCurrentSpecificImpulse( );
    }

    void wrtEngineThrustMagnitude(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::string& engineName )
    {
        std::shared_ptr< system_models::EngineModel > engineModel = engineModelList_.at( engineName );
        partialMatrix( 0, 0 ) += 1.0 /
                ( engineModel->getThrustMagnitudeWrapper( )->getCurrentSpecificImpulse( ) *
                  physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
    }

    virtual bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return false;
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction = std::make_pair( partialFunction, 0 );

        if( parameter->getParameterName( ).first == estimatable_parameters::constant_thrust_magnitude &&
                parameter->getParameterName( ).second.first == body_ &&
                engineModelList_.count( parameter->getParameterName( ).second.second ) != 0 )
        {
            partialFunction = std::make_pair(
                        std::bind( &FromThrustMassRatePartial::wrtEngineThrustMagnitude, this,
                                   std::placeholders::_1,
                                   parameter->getParameterName( ).second.secon ), 1 );

        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::constant_specific_impulse &&
                 parameter->getParameterName( ).second.first == body_ &&
                 engineModelList_.count( parameter->getParameterName( ).second.second ) != 0 )
        {
            partialFunction = std::make_pair(
                        std::bind( &FromThrustMassRatePartial::wrtEngineSpecificImpulse, this,
                                   std::placeholders::_1,
                                   parameter->getParameterName( ).second.second ), 1 );

        }


        return partialFunction;
    }

    virtual void update( const double currentTime )
    {

    }



private:
    std::vector< std::shared_ptr< propulsion::ThrustAcceleration > > thrustAccelerations_;

    std::map< int, std::vector< std::shared_ptr< system_models::EngineModel > > > accelerationBasedThrustSources_;

    std::map< std::string, std::shared_ptr< system_models::EngineModel > > engineModelList_;

};


} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_MASSRATEPARTIAL_H
