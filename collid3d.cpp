#include <FronTier.h>
#include <vector>
#include "collid.h"

//test module for 3d surface
//proximity and collision detection
char *restart_name;
boolean RestartRun;
int RestartStep;
static void propagation_driver(Front*);
static void collision_point_propagate(Front*,POINTER,POINT*,
        			      POINT *newp,HYPER_SURF_ELEMENT *,
        			      HYPER_SURF*,double,double*);
static void initSurfaceState(SURFACE*,const double*);


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
        f_basic.size_of_intfc_state = sizeof(STATE);
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

	const double vel1[] = {-0.1,-0.1,-0.1};
	//initialize velocity function to straight_velocity
	initSurfaceState(surf,vel1);

	center[0] = center[1] = center[2] = 0.20;
            R[0] = R[1] = R[2] = 0.05;
        FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);

	const double vel2[] = {-0.1,-0.1,-0.1};
	initSurfaceState(surf,vel2);

	front.vfunc = NULL;
	//PointPropagationFunction(&front) = fourth_order_point_propagate;
	PointPropagationFunction(&front) = collision_point_propagate;
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
	Tracking_algorithm(front) = STRUCTURE_TRACKING;

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

static void initSurfaceState(
	SURFACE* surf,
	const double* vel)
{
	TRI* tri;
	POINT* p;
	STATE* sl, *sr;
	surf_tri_loop(surf,tri){
	    for (int i = 0; i < 3; ++i){
		p = Point_of_tri(tri)[i];
		sl = (STATE*)left_state(p);
        	sr = (STATE*)right_state(p);
		for (int j = 0; j < 3; ++j){
		    sl->vel[j] = sr->vel[j] = vel[j];
		}
	    }
	}
}

static void collision_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD],s;
        int i, dim = front->rect_grid->dim;

        if (wave_type(oldhs) < MOVABLE_BODY_BOUNDARY)
        {
            for (i = 0; i < dim; ++i)
            {
                Coords(newp)[i] = Coords(oldp)[i];
            }
            return;
        }

	STATE *newsl,*newsr;
        STATE *sl,*sr;
	sl = (STATE*)left_state(oldp);
        sr = (STATE*)right_state(oldp);
        newsl = (STATE*)left_state(newp);
        newsr = (STATE*)right_state(newp);

	for (i = 0; i < dim; ++i)
	    vel[i] = sl->vel[i];

	for (i = 0; i < dim; ++i)
	{
	    newsl->vel[i] = newsr->vel[i] = vel[i];
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
	}
        s = mag_vector(V,dim);
        set_max_front_speed(dim,s,NULL,Coords(newp),front);
}       /* fourth_order_point_propagate */

