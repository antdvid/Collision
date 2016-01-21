#include <FronTier.h>
#include "collid.h"
static void initBalls(Front&);
static void initPlane(Front&);

void initTestModule(Front &front, char* in_name)
{
        FILE* infile = fopen(in_name,"r");
        char mesg[100];
        CursorAfterString(infile,"Enter problem type: ");
        fscanf(infile,"%s",mesg);
        (void) printf("%s\n",mesg);
        if (mesg[0] == 'B' || mesg[0] == 'b')
            initBalls(front);
        else if (mesg[0] == 'P')
	    initPlane(front);
	else
        {
            std::cout << "Unknown problem type" << mesg << std::endl;
            clean_up(ERROR);
        }
}

static void initBalls(Front& front)
{
        double center[MAXD];        // Center of the sphere
        double R[MAXD];
        SURFACE* surf;
        center[0] = center[1] = center[2] = 0.30;
            R[0] = R[1] = R[2] = 0.05;
        FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);

        const double vel1[] = {-0.2,-0.2,-0.2};
        //initialize velocity function to straight_velocity
        initSurfaceState(surf,vel1);

        center[0] = center[1] = center[2] = 0.20;
            R[0] = R[1] = R[2] = 0.05;
        FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,      // refinement level
                        &surf);

        const double vel2[] = {0.2,0.2,0.2};
        initSurfaceState(surf,vel2);

}

static void initPlane(Front& front)
{
        double center[MAXD];        // Center of the sphere
        double R[MAXD];
        SURFACE* surf;
	
	center[0] = center[1] = 0.25;
	center[2] = 0.47;
	double edge[MAXD];
	edge[0] = edge[1] = 0.15;
	edge[2] = 0.01;
	FT_MakeCuboidSurf(&front,center,edge,1,2,
			FIRST_PHYSICS_WAVE_TYPE,
			1,
			&surf);

        const double vel1[] = {0,0,-0.2};
        //initialize velocity function to straight_velocity
        initSurfaceState(surf,vel1);

        center[0] = center[1] = 0.25;
	center[2] = 0.4;
            R[0] = R[1] = R[2] = 0.05;
        FT_MakeEllipticSurf(&front,center,R,
                        1,2,    // negative and positive components
                        MOVABLE_BODY_BOUNDARY,
                        1,      // refinement level
                        &surf);

        const double vel2[] = {0,0,0.2};
        initSurfaceState(surf,vel2);

}


void initSurfaceState(
        SURFACE* surf,
        const double* vel)
{
        TRI* tri;
        POINT* p;
        STATE* sl, *sr;
        surf_tri_loop(surf,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                sl = (STATE*)left_state(p);
                sr = (STATE*)right_state(p);
                for (int j = 0; j < 3; ++j)
                {
                    sl->vel[j] = sr->vel[j] = vel[j];
                    sl->collsnImpulse[j] = sr->collsnImpulse[j] = 0.0;
                    sl->friction[j] = sr->friction[j] = 0.0;
                    sl->collsn_num = sr->collsn_num = 0;
                    sl->x_old[j] = Coords(p)[j];
                }
            }
        }
}

