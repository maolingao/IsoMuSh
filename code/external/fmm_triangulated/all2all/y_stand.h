/********************************************************************/
/* y_stand.h is the declaration part for fast program               */
/********************************************************************/
#ifndef Y_STAND_HDR						/* Make sure that the header file is expanded in */
#define Y_STAND_HDR						/* every .c file at most once.					 */

//#define TrianglesAndEdges				/* When this flag is defined the program create
//											   the files with the data about the triangles
//											   and the edges.							 */

//#define WRITE_PGM						/* When this flag is set the program generates
//											a pgm file after every march.				 */

//#define LEVEL_SET						/* When this flag is defined the program creates
//										   a levelset file.							 */

//#define FULL_OUTPUT						/* This flag is set when we want the full data from
//										   the program, not just the results.			 */

#define NO_DISTANCES					/* Then this flag is defined the program does not
										   generate the distances between the points.	 */

#define NO_MANIFOLD						/* When this flag is defined no maifold files are
										   created.										 */

//#define NO_OLD_U_DATA					/* When this flag is defined the program does not
//										   write the U values of the updating vertexes into
//										   the update files.							 */

//#define NO_NEW_ALIVE_DATA				/* When this flag is defined the program does not create
//										   update records that contain only new alive data with
//										   every other field set to NO_DATA (these record allow
//										   the user of the data to now when every vertex
//										   becomes alive even if it does not update anything
//										   because all its neighbours are alive			 */

#define STOP_ON_SOURCES					/* When this flag is set the program stops the fast
										   march algoritem when the distance of all the sources
										   from the current source is determined.		 */

//#define SMOOTH_GEODESIC					/* When this flag is defined the program calculate
//										   the geodesic lines between the sources, Using the
//										   second order algoritem.						 */
//#define LINEAR_GEODESIC					/* When this flag is defined the program calculate
//										   the geodesic lines between the sources, Using the
//										   first order algoritem.						 */



// #include <windows.h> 
// #include <process.h>


#ifdef __cplusplus
	extern "C++"{
	#include <stdio.h>
	#include <math.h>
	#include <stdlib.h>
	}
#else
	#include <stdio.h>
	#include <math.h>
	#include <stdlib.h>
#endif


#define and 	&&			/* some local defines */
#define or  	||
#define mod 	%
#define not 	!
#define is_not 	!=
#define is	==
#define PI	3.141592654

//#define InSize(N)     (((N)>=0) and ((N)<PicSize))
//#define InFrame(I,J)  (InSize(I) and InSize(J))
#define Sqr(x)  ((x)*(x))
#define Sqrt(x) (sqrt(x))

typedef enum boolean{FALSE = 0, TRUE = 1} boolean;

int	XSIZE,YSIZE;	/* The number of vertexes in the x dimension and the y dimenshion of the grid */
int PICSIZE;			/* The total number of vertixes in the grid. */

			/*big real positive number */
#define MAX_DISTANCE  	9999999    
/*
#define MAX_DISTANCE  	999999999999    
*/

#define MaxTri	 60	/* The maximum number of triangles that a vertex can participate in. */

//#define INT(r) 		((int)((r)+0.5))			/* Round a number to the nearest integer */
//#define IN_RANGEX(x) 	( (x)>=0 && (x)<XSIZE )
//#define IN_RANGEY(y) 	( (y)>=0 && (y)<YSIZE )
//#define IN_RANGE(x,y) 	( IN_RANGEX(x) && IN_RANGEY(y) )
#define ABS(x) 		((x)>0?(x):(-(x)))
#define MAX(x,y) 	((x)>(y)?(x):(y))
#define MIN(x,y) 	((x)>(y)?(y):(x))
#define FAST_MIN(x,y) (((x) > (y)) ? (x) = (y) : 0)	/* Can replace a macro call of the form
													   x = MIN(x, y)					*/
#define MIN4(a,b,c,d) 	MIN(MIN(a,b),MIN(c,d))
//#define MAX4(a,b,c,d) 	MAX(MAX(a,b),MAX(c,d))

#define Alive 0		/* Alive and Far are the values that are placed in the Back Pointer array when the vertex is not */
#define Far   -1	/* on the heap.																					 */


// allocate 2-D matrix: [hsize][vsize] of type _MTYPE
#define matrix_allocate(matrix,_vsize,_hsize,_MTYPE) { \
	int _i; \
	_MTYPE *myTypePtr=NULL; \
	matrix=(_MTYPE **)malloc(_vsize*sizeof(myTypePtr)); \
	if (matrix==NULL) { \
		fprintf(stderr,"Error: Can't allocate memory for 2D array\n"); \
		exit(1); \
	} \
	for (_i=0 ; _i<_vsize; _i++) { \
		matrix[_i]=(_MTYPE *)malloc(_hsize*sizeof(_MTYPE)); \
		if (matrix[_i] ==NULL) { \
			fprintf(stderr,"Error: Can't allocate memory for 2D array\n"); \
			exit(1); \
		} \
	} \
}

// free 2-D matrix: with _vsize lines
#define matrix_free(matrix,_vsize) { \
	int _i; \
	for (_i=0 ; _i<_vsize; _i++) \
		free(matrix[_i]); \
	free(matrix); \
}



struct heap_element { 		/* element of the heap													*/
	double		u;			/* element value (u)													*/
	int 		v; 			/* back pointer to the vertex (the index of the vertex in the array V)	*/
};


#if __cplusplus
	extern "C" int		*BP;			/* back pointer to place in the list */
	extern "C" int		N;				/* number of elements in the heap array */
	extern "C" struct heap_element *a; 	/* heap array */
#else
	extern int			*BP;			/* back pointer to place in the list */
	extern int			N;				/* number of elements in the heap array */
	extern struct heap_element *a; 		/* heap array */
struct heap_element *a; 	/* heap array */
#endif 

struct Triangle { 	/* defines one triangle					*/
	int	Vind[3];	/* Index for the 3 vertices				*/
	double  b[5];   /* surface cooefficients				*/
					/* Du = (2b0x+b2y+b3,2b1y+b2x+b4)		*/
					/* for Vind[0]-Vind[1] = x exis			*/
	boolean Visited;/* Determine whether this vertex was
					   already visited when backtracking to
					   find a smooth geodesic.				*/
	boolean Split;	/* This field insdicates whether on of
					   the edges of this triangle is
					   splitting edge.						*/
};

struct Stencil {		/* virtual numerical connections  */
				/*   true triangles in most cases */
	double	Ctheta;		/* cos(theta) angle between the edges */
	double  Stheta; /* sin(theta)^2 when theta is the angle between the edges */
	int	v1,v2;		/* index of neigboring vertices */
	double	l1,l2;		/* edge length to v1 and to v2 */
	double Ctheta_mul_l1_div_l2_minus_one_mul_2, Ctheta_mul_l2_div_l1_minus_one_mul_2;
								/* Remember the value of the l1 and l2 multiplies
										    by Ctheta minus one	multiplied by 2	*/
	double sqr_l1, sqr_l2; /* The values of l1 and l2 multiplied by themselves. */
						/* Replaces 1 - a2 * (ctheta_mul_mb_div_ma_minus_one_mul_2 + 1) / b2 */ 
	double shortcut1_1, shortcut1_2;
						/* Replaces - Stheta * a2*/ 
	double shortcut2_1, shortcut2_2;

};

struct Vertex { 		/* defines one Numerical Graph vertex */
	double	x,y,z;		/*x,y,z coordinates   */
	double	U;	   	/* U vlaue	*/
	int	Tind[MaxTri];	/* link back to triangles */
					/* MaxTri = Max # triangles at one vertex */
	int	si,vn;		/* number of vertex connections */
					/* si is the #of ST, vn is the index to VN*/
	int ti;			/* The number of triangles this vertex is member of */
	struct  Stencil	ST[MaxTri];	/* numerical connection, 
						updating stenciles  */
	int	VN[3*MaxTri]; /* Neighboring dirrectional vertices indexes*/
				/* Vertex to update */
	int	ST_For_VN[3*MaxTri]; /* The first stenciel of the neighbour with this vertex in it */

	int Split[2];	/* The indexes of the vertexes that split the triangles of that vertex
					   or NO_DATA if there is no vertex.									*/
							 
	boolean source;	  /* Indicates wether the variable is a source in one of the marches
					     (meaning it appears in the xStart, yStart variables).				*/
	boolean current_source;	/* Is it a source that got the U value 0 in the last run of
						 the fast march algoritem.											*/
	int source_num;	  /* In a source vertex this is the index in the xStart and yStart of
					     the coordinates of the vertex.										*/

//	int i,j; // indices to X,Y etc.
#ifdef FULL_OUTPUT
	int u1,u2; // parameter values of the vertex
#endif
//	double E,F,G; // the metric matrix
	double MaxU;
};

#define NoData -1					/* Indicate no data for this field.						*/

#ifdef FULL_OUTPUT

typedef struct split_data{			/* The data about the adding of virtual neighbours for
									   each vertex.											*/
	int virtual_neighbour1_X, virtual_neighbour1_Y;
	int virtual_neighbour2_X, virtual_neighbour2_Y;
	int virtual_neighbour3_X, virtual_neighbour3_Y;
}split_data;

split_data **SplitDataArray;		/* The array where the data about the splitting of the
									   manifold is kept.									*/

typedef struct update_data{			/* One record for the update file.						*/
	double NewAliveX, NewAliveY;	/* The vertex that now becomes alive.					*/
	double UpdatedX, UpdatedY;		/* The vertex that is being updated.					*/
	double ThirdX, ThirdY;			/* The third vertex in the updating triangle or -1.		*/
	double NewU;					/* The new U value of the vertex that is being updated.	*/
#ifndef NO_OLD_U_DATA
	double NewAliveU, ThirdU;		/* The U values of the other vertexes in the stencil	*/
#endif
} update_data;

update_data *UpdateDataArray;		/* Array where all the update data will be kept.		*/

update_data *next_update;			/* Pointer to the next record in the array of the
									   updates.												*/

#endif

struct Point {						/* defines a point in 3D with a U value					*/
	double  x,y,z;					/* x,y,z coordinates IN 3D								*/
	double  U;						/* The U value of the point.							*/
};

struct Vector {						/* defines a vector in 3D								*/
		double  x,y,z;				/* x,y,z coordinates of the vector in 3D.				*/
};

double **U1,			/* The value of the first parametric axes */
	   **U2,			/* The value of the second parametric axes */
	   **X,	**Y,**Z,	/* The values of the coordinates  of every point in the parametric
						   area. */
	   **X1,**X2,		/* The "nigzeret" of X. 1 is by the first parametric axes and 2 is
						   by the other parametric axes	*/
	   **Y1,**Y2,**Z1,**Z2, /* The same as for x, but for y and z */
	   /* E, F, F, G are the values in the metrix. */
	   **E,				/* The sum of the squres of the partial "nigzert" by u1	*/
	   **F,				/* The sum of the multiplications of the partion "nigzert"s*/
	   **G,				/* Like E but for u2. */
	   **COS,			/* The cos of the the axes in the grid point, used to decide wether to
						   split the triangle. */
	   **sqrtE,**sqrtG; /* The squre roots of the values in E and G.*/

double **sources_dis;	/* The distances between every two sources */

#endif /* end of #ifndef Y_STAND_HDR */
