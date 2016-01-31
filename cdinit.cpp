#include <FronTier.h>
#include "collid.h"
static void initBalls(Front&);
static void initPlane(Front&);
static void initString(Front&);
static void initMultipleStrings(Front&);
static CURVE* make3dCurve(Front&,double[][3],int);

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
	else if (mesg[0] == 'S' || mesg[0] == 's')
	    initString(front);
	else if (mesg[0] == 'M')
	    initMultipleStrings(front);
	else
        {
            std::cout << "Unknown problem type: " << mesg << std::endl;
            clean_up(ERROR);
        }
	gview_plot_interface("init",front.interf);
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

static void initString(Front& front){
	double center[MAXD];        // Center of the sphere
        SURFACE* surf;

	//make a cuboid
	center[0] = center[1] = 0.25;
        center[2] = 0.3;
        double edge[MAXD];
        edge[0] = edge[1] = 0.15;
        edge[2] = 0.03;
	FT_MakeCuboidSurf(&front,center,edge,1,2,
                        FIRST_PHYSICS_WAVE_TYPE,
                        1,
                        &surf);	
	const double vel1[] = {0,0,-0.2};
        //initialize velocity function to straight_velocity
        initSurfaceState(surf,vel1);

	//make a string
	double pt[2][3] = {{0.25,0.01,0.25},{0.25,0.49,0.25}};
	CURVE     *curve = make3dCurve(front,pt,STRING_HSBDRY);
}

static void initMultipleStrings(Front& front){
	CURVE *curve;
	//make a string
        double pt1[2][3] = {{0.25,0.01,0.24},{0.25,0.49,0.24}};
	const double vel1[] = {0,0,0.2};
        curve = make3dCurve(front,pt1,STRING_HSBDRY);
	initCurveState(curve,vel1);
	//make a string
        double pt2[2][3] = {{0.01,0.25,0.25},{0.49,0.25,0.25}};
	const double vel2[] = {0,0,-0.2};
        curve = make3dCurve(front,pt2,STRING_HSBDRY);
	initCurveState(curve,vel2);
}

static CURVE* make3dCurve(Front& front,double pt[][3],int hsb_type){
	NODE      *string_nodes[2];
        CURVE     *curve;
        INTERFACE *intfc = front.interf;
        for (int i = 0; i < 2; ++i)
            string_nodes[i] = make_node(Point(pt[i]));
        curve = make_curve(0,0,string_nodes[0],string_nodes[1]);
        hsbdry_type(curve) = hsb_type;
        double spacing = separation(string_nodes[0]->posn,
                                    string_nodes[1]->posn,3);
        double dir[3];
        for (int j = 0; j < 3; ++j)
            dir[j] = (Coords(string_nodes[1]->posn)[j] -
                      Coords(string_nodes[0]->posn)[j])/spacing;
        double *h = computational_grid(intfc)->h;
        int nb = (int)(spacing/(0.25*h[0]));
        spacing /= (double)nb;  
        BOND* b = curve->first; 
        double coords[3];       
        for (int j = 1; j < nb; ++j)
        {
            for (int k = 0; k < 3; ++k)
                coords[k] = Coords(string_nodes[0]->posn)[k] +
                                   j*dir[k]*spacing;
            insert_point_in_bond(Point(coords),b,curve);
            b = b->next;
        }
	return curve;
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

void initCurveState(CURVE* curve, const double* vel)
{
	BOND* b;
	STATE* sl;
	POINT* p[2];
	curve_bond_loop(curve,b)
	{
	    p[0] = b->start;
	    p[1] = b->end;
	    for (int i = 0; i < 2; ++i)
	    {
		sl = (STATE*)left_state(p[i]);
		for (int j = 0; j < 3; ++j)
                {
                    sl->vel[j] = vel[j];
                    sl->collsnImpulse[j] = 0.0;
                    sl->friction[j] = 0.0;
                    sl->collsn_num = 0;
                    sl->x_old[j] = Coords(p[i])[j];
                }	
	    }
	}
}
