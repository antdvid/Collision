#include <FronTier.h>
#include <vector>
#include "collid.h"

//test module for 3d surface
//proximity and collision detection
char *restart_name;
boolean RestartRun;
int RestartStep;
static void propagation_driver(Front*);

int main(int argc, char** argv)
{
	static Front front;
        static F_BASIC_DATA f_basic;
        static LEVEL_FUNC_PACK level_func_pack;
        SURFACE *surf;
	f_basic.dim = 3;
        FT_Init(argc,argv,&f_basic);

        /* Initialize basic computational data */

        f_basic.L[0] = 0.0;     f_basic.L[1] = 0.0;     f_basic.L[2] = 0.0;
        f_basic.U[0] = 0.5;     f_basic.U[1] = 0.5;     f_basic.U[2] = 0.5;
        f_basic.gmax[0] = 50;  f_basic.gmax[1] = 50; f_basic.gmax[2] = 50;
        f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
        f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
        f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
        f_basic.size_of_intfc_state = 0;
	FT_StartUp(&front,&f_basic);
	level_func_pack.pos_component = 2;
	
	restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

	//initialize interface and velocity
        FT_InitIntfc(&front,&level_func_pack);
	double center[MAXD];        // Center of the sphere
        double R[MAXD];
  	center[0] = center[1] = center[2] = 0.30;
            R[0] = R[1] = R[2] = 0.05;	
	FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);
	//initialize velocity function to straight_velocity
	TRANS_PARAMS* trans_params = new TRANS_PARAMS[2];
	trans_params[0].dim = 3;
	trans_params[0].vel[0] = trans_params[0].vel[1] 
			       = trans_params[0].vel[2] = -0.1;
	FT_InitSurfVeloFunc(surf,
			"trans_velocity",
			(POINTER)trans_params,
			translation_vel);

	center[0] = center[1] = center[2] = 0.20;
            R[0] = R[1] = R[2] = 0.05;
        FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);
	//initialize velocity function to zero_velocity
	trans_params[1].dim = 3;
	trans_params[1].vel[0] = trans_params[1].vel[1] 
			       = trans_params[1].vel[2] = 0.0;
	FT_InitSurfVeloFunc(surf,
			"no_velocity",
			trans_params+1,
			translation_vel);
	front.vfunc = NULL;
	PointPropagationFunction(&front) = fourth_order_point_propagate;
	char dname[256];
	sprintf(dname,"%s/intfc",OutName(&front));
	geomview_interface_plot(dname,front.interf,front.rect_grid);

	propagation_driver(&front);
	clean_up(0);
}

static  void propagation_driver(
        Front *front)
{
        double CFL;
	std::vector<std::pair<TRI*,TRI*> > triPairList;
	CollisionSolver *collision_solver = new CollisionSolver3d();

	//collision_solver->setRoundingTolerance(0.1);
	collision_solver->setFabricThickness(0.1);

        front->max_time = 5;
        front->max_step = 100;
        front->print_time_interval = 0.5;
        front->movie_frame_interval = 0.01;

        CFL = Time_step_factor(front) = 0.75;

	Frequency_of_redistribution(front,GENERAL_WAVE) = 100000;
        printf("CFL = %f\n",CFL);
        printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
                Frequency_of_redistribution(front,GENERAL_WAVE));

        if (!RestartRun)
        {
            FT_RedistMesh(front);
            FT_ResetTime(front);

            // Always output the initial interface.
            FT_Save(front);
            FT_Draw(front);

            // This is a virtual propagation to get maximum front 
            // speed to determine the first time step.

            FT_Propagate(front);
            FT_SetTimeStep(front);
            FT_SetOutputCounter(front);
        }
        else
        {
            FT_SetOutputCounter(front);
        }

        FT_TimeControlFilter(front);

        for (;;)
        {
            /* Propagating interface for time step dt */
	    //collision detect and handling
	    triPairList.clear();
	    collision_solver->assembleFromInterface(front->interf);
	    collision_solver->detectProximity();
	    collision_solver->printProximity();
	    char dname[256];
	    sprintf(dname,"%s/intfc",OutName(front));
	    collision_solver->gviewplotPairList(dname);
	    collision_solver->gviewplotPair(dname);
	    geomview_interface_plot(dname,front->interf,front->rect_grid);

            FT_Propagate(front);
            FT_AddTimeStepToCounter(front);

            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

            FT_SetTimeStep(front);

            /* Output section */

            FT_PrintTimeStamp(front);
            if (FT_IsSaveTime(front))
                FT_Save(front);
            if (FT_IsDrawTime(front))
                FT_Draw(front);

            if (FT_TimeLimitReached(front))
                    break;

            FT_TimeControlFilter(front);
        }
	delete collision_solver;
}       /* end propagation_driver */

