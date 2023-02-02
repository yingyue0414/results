/**
 * @file Force.hpp
 * @author Y Ying (yying7@jh.edu)
 * @brief The \code{Energy} defines 7 energy terms and the \code{Force}
 *        class defines corresponding force terms in \code{*gsl_matrix}
 *        type. This file also includes getters and setters for all the
 *        private member variables. Note that all the forces are gsl
 *        matrices with a size of (3,1).
 * @date 2023-01-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <ostream>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <ctime>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "mesh.hpp"
#include "vector_math.hpp" 
#include "parameters.hpp"

/**
 * @brief This class includes 7 components in calculating force
 *        exerted on vertrices: force due to membrane bending,
 *        area constraint, volume constraint, membrane thickness
 *        constraint, lipid tilting, regularization, and spline.
 * 
 */
class Force{
public:
    /****************************************************************/
    /**************************Constructors**************************/
    /****************************************************************/
    /**
     * @brief Construct a new \code{Force} object with all force components set to zero. 
     */
    Force();

    /****************************************************************/
    /*****************************Methods****************************/
    /****************************************************************/
    /**
     * @brief Set \code{forceTotal} to the sum of 7 force components.
     * @return gsl_matrix of calculated total force.
     */
    gsl_matrix* calculateTotalForce();

    /****************************************************************/
    /*****************************Getters****************************/
    /****************************************************************/
    /**
     * @brief Force due to membrane bending.
     * @return gsl_matrix* pointer to the matrix representing forceCurvature
     */
    gsl_matrix* getForceCurvature() const;

    /**
     * @brief Force due to area constraint.
     * @return gsl_matrix* pointer to the matrix representing forceArea
     */
    gsl_matrix* getForceArea() const;

    /**
     * @brief Force due to volume constraint; set 0 if membrane is flat.
     * @return gsl_matrix* pointer to the matrix representing forceVolume
     */
    gsl_matrix* getForceVolume() const;

    /**
     * @brief Force due to membrane thickness change.
     * @return gsl_matrix* pointer to the matrix representing forceThickness
     */
    gsl_matrix* getForceThickness() const;

    /**
     * @brief Force due to lipid tilting.
     * @return gsl_matrix* pointer to the matrix representing forceTilt
     */
    gsl_matrix* getForceTilt() const;

    /**
     * @brief Regularization force; this does not change the
     * geometry of actual limit surface. Instead, it 
     * restrains the shape of each mesh triangle to be
     * as close to equilateral triangle as possible to avoid
     * deformation of the triangular mesh.
     * @return gsl_matrix* pointer to the matrix representing forceRegularization
     */
    gsl_matrix* getForceRegularization() const;

    /**
     * @brief Force due to harmonic bonding between the membrane
     * and the attached spline, which can be a single particle
     * or a lattice. The calculations are defined separately in
     * spline_points.cpp due to it requiring the coordinate
     * information of the spline points.
     * @return gsl_matrix* pointer to the matrix representing forceSpline
     */
    gsl_matrix* getForceHarmonicBond() const;

    /**
     * @brief Total force by summing up the 7 force terms.
     * @return gsl_matrix* pointer to the matrix representing forceTotal
     */
    gsl_matrix* getForceTotal() const;

    /****************************************************************/
    /*****************************Setters****************************/
    /****************************************************************/
    /**
     * @brief Set the force due to membrane bending.
     * 
     * @param forceCurvature A pointer to the gsl_matrix containing the force due to membrane bending.
     */
    void setForceCurvature(gsl_matrix *forceCurvature);

    /**
     * @brief Set the force due to area constraint.
     * 
     * @param forceArea A pointer to the gsl_matrix containing the force due to area constraint.
     */
    void setForceArea(gsl_matrix *forceArea);

    /**
     * @brief Set the force due to volume constraint.
     * 
     * @param forceVolume A pointer to the gsl_matrix containing the force due to volume constraint.
     *                    Should be set to 0 if the membrane is flat.
     */
    void setForceVolume(gsl_matrix *forceVolume);

    /**
     * @brief Set the force due to membrane thickness change.
     * 
     * @param forceThickness A pointer to the gsl_matrix containing the force due to membrane thickness change.
     */
    void setForceThickness(gsl_matrix *forceThickness);

    /**
     * @brief Set the force due to lipid tilting.
     * 
     * @param forceTilt A pointer to the gsl_matrix containing the force due to lipid tilting.
     */
    void setForceTilt(gsl_matrix *forceTilt);

    /**
     * @brief Set the regularization force.
     * 
     * @param forceRegularization A pointer to the gsl_matrix containing the regularization force.
     *                            This does not change the geometry of the actual limit surface.
     *                            Instead, it restrains the shape of each mesh triangle to be
     *                            as close to an equilateral triangle as possible to avoid
     *                            deformation of the triangular mesh.
     */
    void setForceRegularization(gsl_matrix *forceRegularization);

    /**
     * @brief Set the force due to attached spline.
     * 
     * @param forceSpline A pointer to the gsl_matrix containing the force due to attached spline.
     *                    The spline can be a single particle or a lattice. 
     *                    The calculations are defined separately in spline_points.cpp 
     *                    due to it requiring the coordinate information of the spline points.
     */
    void setForceHarmonicBond(gsl_matrix *forceSpline);

    /**
     * @brief Set the total force.
     * 
     * @param forceTotal A pointer to the gsl_matrix containing the total force. 
     *                   It is obtained by summing up the 7 force terms.
     */
    void setForceTotal(gsl_matrix *forceTotal);
    
private:
    
    gsl_matrix* _forceCurvature; // Force due to membrane bending.
    gsl_matrix* _forceArea; // Force due to area constraint.
    gsl_matrix* _forceVolume; // Force due to volume constraint; set 0 if membrane is flat.
    gsl_matrix* _forceThickness; // Force due to membrane thickness change.
    gsl_matrix* _forceTilt; // Force due to lipid tilting.
    gsl_matrix* _forceRegularization; // Regularization force; this does not change the
                                      // geometry of actual limit surface. Instead, it 
                                      // restrains the shape of each mesh triangle to be
                                      // as close to equilateral triangle as possible to avoid
                                      // deformation of the triangular mesh.
    gsl_matrix* _forceSpline; // Force due to attached spline, which can be a single particle
                              // or a lattice. The calculations are defined separately in
                              // spline_points.cpp due to it requiring the coordinate
                              // information of the spline points.
    gsl_matrix* _forceTotal; // Total force by summing up the 7 force terms.

}

/****************************************************************/
/*****************************ostream****************************/
/****************************************************************/
std::ostream& operator<<(std::ostream &stream, const Force &force);