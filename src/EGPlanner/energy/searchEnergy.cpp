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

#define HAND_POSE_PITCH_FILTER_UPPER                         (1.0)          /* 60 degrees cone up and down */
#define HAND_POSE_PITCH_FILTER_LOWER                         (-1.0)

#define HAND_POSE_ROLL_FILTER_UPPER                         (1.0)          /* 60 degrees cone up and down */
#define HAND_POSE_ROLL_FILTER_LOWER                         (-1.0)

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

    double hand_roll , hand_pitch , hand_yaw ;

    mat3 hand_rotation_matrix ;
    hand_rotation.ToRotationMatrix(hand_rotation_matrix);
    hand_rotation_matrix.ToEulerAngles(hand_roll,hand_pitch,hand_yaw);

    Quaternion object_rotation = objTran.rotation();
    vec3 object_translation = objTran.translation();

    bool grasp_out_of_limit = false ;
    bool grasp_x_axis_exceeded = false ;
    bool grasp_pitch_exceeded = false ;
    bool grasp_roll_exceeded = false ;
    double position_violation_penalty = 0 , pitch_violation_penalty = 0 , roll_violation_penalty = 0 ;

    double hand_pitch_upper_limit , hand_pitch_lower_limit , hand_roll_upper_limit , hand_roll_lower_limit ;


    if(hand_yaw < 0)
    {
        /* If Yaw is decreasing from 0 to -3.14 then Pitch and Roll will increase being Pitch increasing at half 
           the rate of Yaw and Roll increasing proportionately */
        hand_pitch_upper_limit = HAND_POSE_PITCH_FILTER_UPPER - (hand_yaw / 2) ;
        hand_pitch_lower_limit = HAND_POSE_PITCH_FILTER_LOWER + (hand_yaw / 2) ;
        hand_roll_upper_limit = HAND_POSE_ROLL_FILTER_UPPER - hand_yaw ;
        hand_roll_lower_limit = HAND_POSE_ROLL_FILTER_LOWER + hand_yaw ;
    }
    else
    {
        /* If Yaw is increasing from -3.14 to 0 then Pitch and Roll will decrease being Pitch decreasing at half 
           the rate of Yaw and Roll decreasing proportionately */
        hand_pitch_upper_limit = HAND_POSE_PITCH_FILTER_UPPER - (hand_yaw / 2) ;
        hand_pitch_lower_limit = HAND_POSE_PITCH_FILTER_LOWER + (hand_yaw / 2) ;
        hand_roll_upper_limit = HAND_POSE_ROLL_FILTER_UPPER - hand_yaw ;
        hand_roll_lower_limit = HAND_POSE_ROLL_FILTER_LOWER + hand_yaw ;
    }


    /*hand_pitch_upper_limit = HAND_POSE_PITCH_FILTER_UPPER + hand_yaw ;
    hand_pitch_lower_limit = HAND_POSE_PITCH_FILTER_UPPER - hand_yaw ;
    hand_roll_upper_limit = HAND_POSE_ROLL_FILTER_UPPER + hand_yaw ;
    hand_roll_lower_limit = HAND_POSE_ROLL_FILTER_LOWER - hand_yaw ;*/

    if(hand_translation.x() > object_translation.x())
    {
        grasp_out_of_limit = true ;
        grasp_x_axis_exceeded = true ;
        position_violation_penalty = hand_translation.x() - object_translation.x() ;
    }

    if((hand_pitch >= HAND_POSE_PITCH_FILTER_LOWER) && (hand_pitch <= HAND_POSE_PITCH_FILTER_UPPER))
    {
        grasp_out_of_limit = true ;
        grasp_pitch_exceeded = true ;
        if(hand_pitch >= 0)
        {
            pitch_violation_penalty = ((hand_pitch_upper_limit*2) - hand_pitch) ;
        }
        else
        {
            pitch_violation_penalty = ((hand_pitch_upper_limit*2) + hand_pitch) ;
        }
    }

    if(!((hand_roll >= HAND_POSE_ROLL_FILTER_LOWER) && (hand_roll <= HAND_POSE_ROLL_FILTER_UPPER)))
    {
        grasp_out_of_limit = true ;
        grasp_roll_exceeded = true ;
        if(hand_roll >= 0)
        {
            roll_violation_penalty = ((hand_roll_upper_limit*2) - hand_roll) ;
        }
        else
        {
            roll_violation_penalty = ((hand_roll_upper_limit*2) + hand_roll) ;
        }
    }

    if ( !state->execute() || !legal() ) {
        isLegal = false;
        stateEnergy = 0;
    } else {
        isLegal = true;
        if(grasp_out_of_limit)
        {
            if(grasp_x_axis_exceeded && grasp_pitch_exceeded && grasp_roll_exceeded) {
                stateEnergy = energy() + position_violation_penalty + pitch_violation_penalty + roll_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Position + Pitch + Roll Penalty error " << position_violation_penalty + pitch_violation_penalty + roll_violation_penalty );
            }
            else if(grasp_x_axis_exceeded && grasp_pitch_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + pitch_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Position + Pitch Penalty Penalty error " << position_violation_penalty + pitch_violation_penalty);
            }
            else if(grasp_pitch_exceeded && grasp_roll_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty + roll_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Pitch + Roll Penalty error " << pitch_violation_penalty + roll_violation_penalty);
            }
            else if(grasp_x_axis_exceeded && grasp_roll_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty + roll_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Position + Roll Penalty error " << position_violation_penalty + roll_violation_penalty);
            }
            else if(grasp_x_axis_exceeded)
            {
                stateEnergy = energy() + position_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Position Penalty error " << position_violation_penalty);
            }
            else if(grasp_pitch_exceeded)
            {
                stateEnergy = energy() + pitch_violation_penalty  ;   //Adding Penalty and increasing energy
                DBGA("Pitch Penalty error " << pitch_violation_penalty);
            }
            else if(grasp_roll_exceeded)
            {
                stateEnergy = energy() + roll_violation_penalty ;   //Adding Penalty and increasing energy
                DBGA("Pitch + Roll Penalty error " << roll_violation_penalty);
            }
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




