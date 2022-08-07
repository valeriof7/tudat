/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      A. Dobrovolskis (1996), "Inertia of Any Polyhedron", Icarus, 124 (243), 698-704
 */

#ifndef TUDATBUNDLE_POLYHEDRONFUNTIONS_H
#define TUDATBUNDLE_POLYHEDRONFUNTIONS_H

#include <Eigen/Core>

namespace tudat
{
namespace polyhedron_utilities
{

/*! Checks if the provided polyhedron settings are valid.
 *
 * Checks if the provided polyhedron settings are valid. It verifies that the provided matrices with the vertices coordinates
 * and the vertices defining each facet have valid dimensions. Throwns an error if invalid dimensions.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet).
 */
void checkValidityOfPolyhedronSettings( const Eigen::MatrixXd& verticesCoordinates,
                                        const Eigen::MatrixXi& verticesDefiningEachFacet);

/*! Computes the volume of a polyhedron.
 *
 * Computes the volume of a polyhedron, according to Dobrovolskis (1996), section 3.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet).
 * @return Volume.
 */
double computeVolume ( const Eigen::MatrixXd& verticesCoordinates,
                       const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Computes the centroid a polyhedron.
 *
 * Computes the volume of a polyhedron, according to Dobrovolskis (1996), section 4. When using the polyhedron as a
 * constant density gravity model, the centroid coincides with the center of mass.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet).
 * @return Position of centroid.
 */
Eigen::Vector3d computeCentroidPosition (const Eigen::MatrixXd& verticesCoordinates,
                                         const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Adjusts the centroid of the polyhedron to the desired value.
 *
 * Adjusts the centroid of the polyhedron to the desired value. When using the polyhedron as a constant density gravity
 * model, the centroid coincides with the center of mass; this function might be useful e.g. to ensure that the center
 * of mass coincides with the origin of some body fixed frame.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet).
 * @param desiredCentroid Desired position of the centroid wrt the frame in which the vertices were defined.
 * @return Corrected coordinates of vertices.
 */
Eigen::MatrixXd modifyCentroidPosition ( Eigen::MatrixXd verticesCoordinates,
                                         const Eigen::MatrixXi& verticesDefiningEachFacet,
                                         const Eigen::Vector3d desiredCentroid );

} // namespace polyhedron_utilities
} // namespace tudat

#endif //TUDATBUNDLE_POLYHEDRONFUNTIONS_H
