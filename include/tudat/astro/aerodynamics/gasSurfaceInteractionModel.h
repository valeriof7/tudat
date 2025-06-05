/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
#ifndef TUDAT_GASSURFACEINTERACTIONMODEL_H
#define TUDAT_GASSURFACEINTERACTIONMODEL_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{

class GasSurfaceInteractionModel
{
public:
GasSurfaceInteractionModel( const GasSurfaceInteractionModelType modelType,
                            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                            const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction,
                            const int maximumNumberOfPixels, 
                            const std::function< Eigen::Vector3d( ) > airSpeedVectorFunction, 
                            const std::function< double( ) > freeStreamTemperatureFunction ):
                            modelType_( modelType), allPanels_( allPanels ), 
                            rotationToBodyFrameFunction_( rotationToBodyFrameFunction ),
                            maximumNumberOfPixels_( maximumNumberOfPixels ),
                            airSpeedVectorFunction_( airSpeedVectorFunction ),
                            freeStreamTemperatureFunction_( freeStreamTemperatureFunction )                            
{ 
    totalNumberOfPanels_ = allPanels_.size( );
    surfacePanelCosines_.resize( totalNumberOfPanels_ );
    if ( maximumNumberOfPixels_ > 2 )
    {
        aerodynamicSelfShadowing_ = std::make_shared< tudat::system_models::SelfShadowing >( allPanels, maximumNumberOfPixels_ );
    }
    else
    {
        aerodynamicSelfShadowing_ = nullptr;
    }
    unityIlluminationFraction_ = std::vector< double >( totalNumberOfPanels_, 1.0);
}

virtual Eigen::Vector3d computeAerodynamicCoefficients( ) = 0;

virtual ~GasSurfaceInteractionModel( ) = default;

void updateMembers( );

double getReferenceArea( ) const
{
    return referenceArea_;
}

std::vector< double >& getSurfacePanelCosines( )
{
    return surfacePanelCosines_;
}

std::vector< double >& getIlluminatedPanelFractions( )
{
    return illuminatedPanelFractions_;
}

protected:
GasSurfaceInteractionModelType modelType_;

const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels_;

const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction_;

std::shared_ptr< tudat::system_models::SelfShadowing > aerodynamicSelfShadowing_;

int maximumNumberOfPixels_;

double referenceArea_;

std::function< Eigen::Vector3d( ) > airSpeedVectorFunction_;

std::function< double( ) > freeStreamTemperatureFunction_;

Eigen::Vector3d incomingDirection_;

std::vector< double > unityIlluminationFraction_;

std::vector< double > illuminatedPanelFractions_;

std::vector< double > surfacePanelCosines_;

int totalNumberOfPanels_;

double speedRatio_;

};

class NewtonGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
NewtonGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                  const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction,
                                  const int maximumNumberOfPixels, 
                                  const std::function< Eigen::Vector3d( ) > airSpeedVectorFunction, 
                                  const std::function< double( ) > freeStreamTemperatureFunction ):
                                  GasSurfaceInteractionModel( newton, allPanels, rotationToBodyFrameFunction,
                                    maximumNumberOfPixels, airSpeedVectorFunction, freeStreamTemperatureFunction )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class StorchGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
StorchGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                  const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction,
                                  const int maximumNumberOfPixels, 
                                  const std::function< Eigen::Vector3d( ) > airSpeedVectorFunction, 
                                  const std::function< double( ) > freeStreamTemperatureFunction ):
                                  GasSurfaceInteractionModel( storch, allPanels, rotationToBodyFrameFunction,
                                    maximumNumberOfPixels, airSpeedVectorFunction, freeStreamTemperatureFunction )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class SentmanGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
SentmanGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                   const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction,
                                   const int maximumNumberOfPixels, 
                                   const std::function< Eigen::Vector3d( ) > airSpeedVectorFunction, 
                                   const std::function< double( ) > freeStreamTemperatureFunction ):
                                   GasSurfaceInteractionModel( sentman, allPanels, rotationToBodyFrameFunction,
                                     maximumNumberOfPixels, airSpeedVectorFunction, freeStreamTemperatureFunction )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class CookGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
CookGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                const std::function< Eigen::Quaterniond( ) > rotationToBodyFrameFunction,
                                const int maximumNumberOfPixels, 
                                const std::function< Eigen::Vector3d( ) > airSpeedVectorFunction, 
                                const std::function< double( ) > freeStreamTemperatureFunction ):
                                GasSurfaceInteractionModel( cook, allPanels, rotationToBodyFrameFunction,
                                  maximumNumberOfPixels, airSpeedVectorFunction, freeStreamTemperatureFunction )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_GASSURFACEINTERACTIONMODEL_H