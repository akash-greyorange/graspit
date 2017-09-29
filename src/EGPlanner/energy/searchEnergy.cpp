//######################################################################
//
// GraspIt!
// Copyright (C) 2002-2009  Columbia University in the City of New York.
// All rights reserved.
//
// GraspIt! is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GraspIt! is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GraspIt!.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s): Matei T. Ciocarlie
//
// $Id: searchEnergy.cpp,v 1.42 2009/09/13 19:57:38 hao Exp $
//
//######################################################################

#include "searchEnergy.h"

#include <time.h>

#include "robot.h"
#include "barrett.h"
#include "body.h"
#include "grasp.h"
#include "contact/contact.h"
#include "world.h"
#include "quality.h"
#include "searchState.h"
#include "graspitCore.h"
#include "ivmgr.h"
#include "matrix.h"

#include "EGPlanner/energy/contactEnergy.h"
#include "EGPlanner/energy/potentialQualityEnergy.h"
#include "EGPlanner/energy/guidedPotentialQualityEnergy.h"
#include "EGPlanner/energy/autograspQualityEnergy.h"
#include "EGPlanner/energy/guidedAutoGraspEnergy.h"
#include "EGPlanner/energy/dynamicAutograspEnergy.h"
#include "EGPlanner/energy/compliantEnergy.h"
#include "EGPlanner/energy/strictAutoGraspEnergy.h"

//#define GRASPITDBG
#include "debug.h"

//#define PROF_ENABLED
#include "profiling.h"

#include <matvec3D.h>

PROF_DECLARE(QS);

#define HAND_POSE_ROLL_RANGE_1_UPPER                            (0.0)
#define HAND_POSE_ROLL_RANGE_1_LOWER                            (-0.30)

#define HAND_POSE_ROLL_RANGE_2_UPPER                            (0.30)
#define HAND_POSE_ROLL_RANGE_2_LOWER                            (0.0)

#define HAND_POSE_ROLL_RANGE_3_UPPER                            (-2.84)
#define HAND_POSE_ROLL_RANGE_3_LOWER                            (-3.14)

#define HAND_POSE_ROLL_RANGE_4_UPPER                            (3.14)
#define HAND_POSE_ROLL_RANGE_4_LOWER                            (2.84)


#define HAND_POSE_PITCH_RANGE_1_UPPER                            (1.57)          
#define HAND_POSE_PITCH_RANGE_1_LOWER                            (1.27)

#define HAND_POSE_PITCH_RANGE_2_UPPER                           (1.87)
#define HAND_POSE_PITCH_RANGE_2_LOWER                           (1.57)

#define HAND_POSE_PITCH_RANGE_3_UPPER                           (-1.27)   
#define HAND_POSE_PITCH_RANGE_3_LOWER                           (-1.57)

#define HAND_POSE_PITCH_RANGE_4_UPPER                           (-1.57)
#define HAND_POSE_PITCH_RANGE_4_LOWER                           (-1.87)


#define HAND_POSE_YAW_RANGE_1_UPPER                             (0.0)
#define HAND_POSE_YAW_RANGE_1_LOWER                             (-0.30)

#define HAND_POSE_YAW_RANGE_2_UPPER                             (0.30)    
#define HAND_POSE_YAW_RANGE_2_LOWER                             (0.0)

#define HAND_POSE_YAW_RANGE_3_UPPER                             (3.14)
#define HAND_POSE_YAW_RANGE_3_LOWER                             (2.84)

#define HAND_POSE_YAW_RANGE_4_UPPER                             (-2.84)
#define HAND_POSE_YAW_RANGE_4_LOWER                             (-3.14)

#define OBJECT_GRASPING_DISTANCE                                (180)


//todo move this out of here
const double unbalancedForceThreshold = 1.0e10;

SearchEnergy::SearchEnergy()
{
    mHand = NULL;
    mObject = NULL;
    mType = ENERGY_CONTACT; //default
    mContactType = CONTACT_LIVE; //default
    mVolQual = NULL;
    mEpsQual = NULL;
    mDisableRendering = true;
    mOut = NULL;
    mThreshold = 0;
    mAvoidList = NULL;
}

void
SearchEnergy::createQualityMeasures()
{
    if (mVolQual) delete mVolQual;
    if (mEpsQual) delete mEpsQual;
    mVolQual = new QualVolume( mHand->getGrasp(), QString("SimAnn_qvol"),"L1 Norm");
    mEpsQual = new QualEpsilon( mHand->getGrasp(), QString("SimAnn_qeps"), "L1 Norm");
    DBGP("Qual measures created");
}

void
SearchEnergy::setHandAndObject(Hand *h, Body *o)
{
    if (mHand != h) {
        mHand = h;
        createQualityMeasures();
    }
    mObject = o;
}

SearchEnergy::~SearchEnergy()
{
    delete mVolQual;
    delete mEpsQual;
}

bool
SearchEnergy::legal() const
{	
    //hack for iros09
    //the compliant planners do their own checks
    if (mType == ENERGY_COMPLIANT || mType == ENERGY_DYNAMIC) return true;

    //no check at all
    //return true;

    //full collision detection
    //if the hand is passed as an argument, this should only check for collisions that
    //actually involve the hand
    bool collision = mHand->getWorld()->noCollision(mHand) ;
    /*if(collision)
    {
        DBGA("Hand is not in collision with environment");
    }
    else
    {
        DBGA("Hand is in collision with environment");
    }*/

    return collision ;

    /*
    //check only palm
    if ( mHand->getWorld()->getDist( mHand->getPalm(), mObject) <= 0) return false;
    return true;
    */
}

void
SearchEnergy::analyzeCurrentPosture(Hand *h, Body *o, bool &isLegal, double &stateEnergy, bool noChange)
{
    setHandAndObject(h,o);

    if (noChange) {
        h->saveState();
    }

    if ( !legal() ) {
        isLegal = false;
        stateEnergy = 0;
    } else {
        isLegal = true;
        stateEnergy = energy();
    }

    if (noChange) {
        h->restoreState();
    }
}

void SearchEnergy::analyzeState(bool &isLegal, double &stateEnergy, const GraspPlanningState *state, bool noChange)
{
    if (mAvoidList) {
        std::list<GraspPlanningState*>::const_iterator it; int i=0;
        for (it = mAvoidList->begin(); it!=mAvoidList->end(); it++){
            if ( (*it)->distance(state) < mThreshold ) {
                isLegal = false;
                stateEnergy = 0.0;
                DBGP("State rejected; close to state " << i);
                return;
            }
            i++;
        }
    }

    Hand *h = state->getHand();
    setHandAndObject( h, state->getObject() );
    h->saveState();
    transf objTran = state->getObject()->getTran();

    bool render = h->getRenderGeometry();
    if( mDisableRendering) {
        h->setRenderGeometry(false);
    }

    transf hand_position = state->getTotalTran();
    Quaternion hand_rotation = hand_position.rotation();
    vec3 hand_translation = hand_position.translation();

    /*Quaternion sample_rot(-0.463032 ,0.542077,0.472774,0.517918);
    mat3 hand_rotation_matrix_sample ;
    sample_rot.ToRotationMatrix(hand_rotation_matrix_sample);
    double hand_roll_sample , hand_pitch_sample , hand_yaw_sample ;
    hand_rotation_matrix_sample.ToEulerAngles(hand_roll_sample,hand_pitch_sample,hand_yaw_sample);
    DBGA("Fixed Hand Roll:" << hand_roll_sample << "Fixed Hand Pitch:" << hand_pitch_sample << "Fixed Hand Yaw:" << hand_yaw_sample);*/

    double hand_roll , hand_pitch , hand_yaw ;

    mat3 hand_rotation_matrix ;
    hand_rotation.ToRotationMatrix(hand_rotation_matrix);
    hand_rotation_matrix.ToEulerAngles(hand_yaw,hand_pitch,hand_roll);
    if(hand_roll > 0)
    {
        hand_roll = 3.14 - hand_roll ;
    }
    else
    {
        hand_roll = -3.14 + hand_roll ;
    }

    Quaternion object_rotation = objTran.rotation();
    vec3 object_translation = objTran.translation();

    bool grasp_out_of_limit = false ;
    bool grasp_x_axis_exceeded = false ;
    bool grasp_pitch_exceeded = false ;
    bool grasp_roll_exceeded = false ;
    bool grasp_yaw_exceeded = false ;
    double position_violation_penalty = 0 , pitch_violation_penalty = 0 , roll_violation_penalty = 0 , yaw_violation_penalty = 0 ;
    
    if(hand_translation.x() > (object_translation.x() - OBJECT_GRASPING_DISTANCE))
    {
        grasp_out_of_limit = true ;
        grasp_x_axis_exceeded = true ;
        position_violation_penalty = (hand_translation.x() - (object_translation.x() - OBJECT_GRASPING_DISTANCE)) * 100 ;
    }
    /*else 
    {
        if(!(((hand_pitch >= HAND_POSE_PITCH_RANGE_1_LOWER) && (hand_pitch <= HAND_POSE_PITCH_RANGE_1_UPPER)) || ((hand_pitch >= HAND_POSE_PITCH_RANGE_2_LOWER) &&
            (hand_pitch <= HAND_POSE_PITCH_RANGE_2_UPPER)) || ((hand_pitch >= HAND_POSE_PITCH_RANGE_3_LOWER) && (hand_pitch <= HAND_POSE_PITCH_RANGE_3_UPPER))
            || ((hand_pitch >= HAND_POSE_PITCH_RANGE_4_LOWER) && (hand_pitch <= HAND_POSE_PITCH_RANGE_4_UPPER))))
        {
            grasp_out_of_limit = true ;
            grasp_pitch_exceeded = true ;
            if(hand_pitch >= 0)
            {
                pitch_violation_penalty = (1.57 - hand_pitch) * 100 ;
            }
            else
            {
                pitch_violation_penalty = (1.57 + hand_pitch) * 100 ;
            }
        }
        else
        {
            if(!(((hand_yaw >= HAND_POSE_YAW_RANGE_1_LOWER) && (hand_yaw <= HAND_POSE_YAW_RANGE_1_UPPER)) || ((hand_yaw >= HAND_POSE_YAW_RANGE_2_LOWER) &&
            (hand_yaw <= HAND_POSE_YAW_RANGE_2_UPPER)) || ((hand_yaw >= HAND_POSE_YAW_RANGE_3_LOWER) && (hand_yaw <= HAND_POSE_YAW_RANGE_3_UPPER)) ||
            ((hand_yaw >= HAND_POSE_YAW_RANGE_4_LOWER) && (hand_yaw <= HAND_POSE_YAW_RANGE_4_UPPER))))
            {
                grasp_out_of_limit = true ;
                grasp_yaw_exceeded = true ;
                if(hand_yaw >= 0)
                {
                    yaw_violation_penalty = hand_yaw * 100 ;
                }
                else
                {
                    yaw_violation_penalty -= hand_yaw * 100 ;
                }
            }

            if(!(((hand_roll >= HAND_POSE_ROLL_RANGE_1_LOWER) && (hand_roll <= HAND_POSE_ROLL_RANGE_1_UPPER)) || ((hand_roll >= HAND_POSE_ROLL_RANGE_2_LOWER) && 
                (hand_roll <= HAND_POSE_ROLL_RANGE_2_UPPER)) || ((hand_roll >= HAND_POSE_ROLL_RANGE_3_LOWER) && (hand_roll <= HAND_POSE_ROLL_RANGE_3_UPPER))
                || ((hand_roll >= HAND_POSE_ROLL_RANGE_4_LOWER) && (hand_roll <= HAND_POSE_ROLL_RANGE_4_UPPER))))
            {
                grasp_out_of_limit = true ;
                grasp_roll_exceeded = true ;
                if(hand_roll >= 0)
                {
                    roll_violation_penalty = (3.14 - hand_roll) * 100 ;
                }
                else
                {
                    roll_violation_penalty -= ( -3.14 + hand_roll) * 100 ;
                }
            }
        }
    }*/

    if ( !state->execute() || !legal() ) {
        isLegal = false;
        stateEnergy = 0;
    } else {
        isLegal = true;
        if(grasp_out_of_limit)
        {
            /*if(grasp_x_axis_exceeded && grasp_pitch_exceeded && grasp_roll_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + pitch_violation_penalty + roll_violation_penalty + yaw_violation_penalty ;
            }
            else if(grasp_x_axis_exceeded && grasp_pitch_exceeded && grasp_roll_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + pitch_violation_penalty + roll_violation_penalty ;
            }
            else if(grasp_x_axis_exceeded && grasp_roll_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + roll_violation_penalty + yaw_violation_penalty ;
            }
            else if(grasp_x_axis_exceeded && grasp_pitch_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + pitch_violation_penalty + yaw_violation_penalty ;
            }
            else if(grasp_pitch_exceeded && grasp_roll_exceeded && grasp_yaw_exceeded) 
            {
                stateEnergy = energy() + pitch_violation_penalty + roll_violation_penalty + yaw_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_x_axis_exceeded && grasp_pitch_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty + position_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_pitch_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty + yaw_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_pitch_exceeded && grasp_roll_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty + roll_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_roll_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + roll_violation_penalty + yaw_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_x_axis_exceeded && grasp_roll_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + roll_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_x_axis_exceeded && grasp_yaw_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + yaw_violation_penalty ;   //Adding Penalty and increasing energy
            }*/
            if(grasp_x_axis_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty ;   //Adding Penalty and increasing energy
            }
            /*else if(grasp_pitch_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty  ;   //Adding Penalty and increasing energy
            }
            else if(grasp_roll_exceeded)
            {
                stateEnergy = energy() + roll_violation_penalty ;   //Adding Penalty and increasing energy
            }
            else if(grasp_yaw_exceeded)
            {
                stateEnergy = energy() + yaw_violation_penalty ;    //Adding Penalty and increasing energy
            }*/
        }
        else
        {
            stateEnergy = energy();
        }
        //isLegal = true;
        //stateEnergy = energy();
    }

    if (noChange || !isLegal) {
        h->restoreState();
        state->getObject()->setTran(objTran);
    }

    if (render && mDisableRendering) h->setRenderGeometry(true);
    return;
}


double SearchEnergy::getEpsQual(){
    mHand->getWorld()->findAllContacts();
    mHand->getWorld()->updateGrasps();
    return mEpsQual->evaluate();
}

double SearchEnergy::getVolQual(){
    mHand->getWorld()->findAllContacts();
    mHand->getWorld()->updateGrasps();
    return mVolQual->evaluate();
}

SearchEnergy * SearchEnergy::getSearchEnergy(SearchEnergyType type)
{
    SearchEnergy *se;

    switch (type)
    {
    case ENERGY_CONTACT:
        se = new ContactEnergy();
        break;
    case ENERGY_POTENTIAL_QUALITY:
        se =  new PotentialQualityEnergy();
        break;
    case ENERGY_AUTOGRASP_QUALITY:
        se =  new AutoGraspQualityEnergy();
        break;
    case ENERGY_CONTACT_QUALITY:
        se =  new GuidedPotentialQualityEnergy();
        break;
    case ENERGY_GUIDED_AUTOGRASP:
        se =  new GuidedAutoGraspQualityEnergy();
        break;
    case ENERGY_STRICT_AUTOGRASP:
        se =  new StrictAutoGraspEnergy();
        break;
    case ENERGY_COMPLIANT:
        se =  new CompliantEnergy();
        break;
    case ENERGY_DYNAMIC:
        se =  new DynamicAutoGraspEnergy();
        break;
    default:
        std::cout << "INVALID SEARCH ENERGY TYPE: " <<  type << std::endl;
        return NULL;
    }

    se->setType(type);
    return se;
}




