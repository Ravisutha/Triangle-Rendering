/************************************************************ 	
 *	Author 		: Ravisutha Sakrepatna Srinivasamurthy
 *	Assignment	: Image Rendering 	
 *************************************************************/
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Macros */
#define ERR -1
#define SUCCESS 1
#define COLS 256
#define ROWS 256

typedef struct vector
{
		float x, y, z;
}vec;

/* Function prototype */
int get_num_ver_faces (int *faces, int *vertices, FILE *fp);
int get_co_ordinates (int vertices, int faces, float *x, float *y, float *z, int *a, int *b, int *c, FILE *fp);
int print_faces (int faces, int *a, int *b, int *c);
int print_vertices (int vertices, float *x, float *y, float *z);
int find_min_max (int vertices, float *x, float *y, float *z, vec *min, vec *max);
int get_center (vec min, vec max, vec *center);
float get_extent (vec min, vec max);
void initialize_vector (vec *camera, int x, int y, int z);
void rotation(int x, float ***r_x, float theta); 
void vector_rotation (vec *camera, vec *up, float **A);
void print_matrix (float **r);
float covert_rad (float theta);
void move_scale (float extent, vec *camera, vec *center);
void cross_product (vec A, vec B, vec *C);
void mul (float num, vec *A);
void sum (vec A, vec B, vec *C);
void sub (vec A, vec B, vec *C);
void dot_product (vec A, vec B, float *C);
void determine_3D (vec *left, vec *right, vec *up, vec *bottom, vec *top, vec *topleft, vec *center, vec camera, float e);
float norm (vec A);
void print_vector (vec A);
void print_cord (vec left, vec right, vec up, vec bottom, vec topleft, vec center, vec camera);
int matrix_multiplication (float **A, float **B, float **C, int num);
int vector_matrix_multiplication (float *V, float **A, int index, float *R);
void vector_coordinates (vec *image, vec right, vec left, vec bottom, vec top, vec topleft, int cols, int rows);
void plane_equation (vec *vert, int *a, int *b, int *c, int faces, vec *ABC, float *D);
void write_ppm (vec *intersect, int faces, int cols, int rows, int *skip, int num, vec *vert, int *a, int *b, int *c, FILE *fp);
void intersection_point (float *dot1, float *dot2, float *dot3, vec intersect, vec v0, vec v1, vec v2);
int find_intersect (int vertices, float *n, float *d, int *skip, vec *image, vec camera, vec *intersect);
void find_distance (float *n, float *d, int faces, vec *ABC, float *D, vec camera, vec *image);
void find_pixel (vec *vert, vec camera, int *a, int *b, int *c, int rows, int cols, int faces, FILE *fp, vec left, vec right, vec top, vec bottom, vec topleft);

int main (int argc, char **argv)
{
		/* Declarations */
		FILE *fp, *out_fp;;
		int faces, vertices, i;
		float *x, *y, *z;
		float theta1, theta2, theta3;
		int *a, *b, *c, cols = COLS, rows = ROWS;
		vec min, max, center, *vert;
		vec camera, up, left, right, top, bottom, topleft, *intersect;
		float **r_x, **r_y, **r_z;
		float e, **A, **B, **C, *D, *n, *d;

		if (argc != 5)
		{
				printf ("Input the file name\n");
				return ERR;
		}

		if ((out_fp = fopen ("rendered.ppm", "wb")) == 0)
		{
				perror ("Output:");
				return ERR;
		}

		/* Memory allocation */
		A = calloc (3, sizeof (float *));
		B = calloc (3, sizeof (float *));
		C = calloc (3, sizeof (float *));

		for (i = 0; i < 3; i++)
		{
				A[i] = calloc (3, sizeof (float));
				B[i] = calloc (3, sizeof (float));
				C[i] = calloc (3, sizeof (float));
		}

		/* Get theta in degrees */
		theta1 = atof (argv[2]);
		theta1 = covert_rad (theta1);
		theta2 = atof (argv[3]);
		theta2 = covert_rad (theta2);
		theta3 = atof (argv[4]);
		theta3 = covert_rad (theta3);

		if ((fp = fopen (argv[1], "r")) == NULL)
		{
				perror ("File opening:");
				return ERR;
		}

		/* Get number of vertices and faces */
		get_num_ver_faces (&faces, &vertices, fp);

		/* Allocate memory */
		x = calloc (vertices, sizeof (float));
		y = calloc (vertices, sizeof (float));
		z = calloc (vertices, sizeof (float));
		a = calloc (faces, sizeof (int));
		b = calloc (faces, sizeof (int));
		c = calloc (faces, sizeof (int));

		/* Get the (x, y, z) of the vertices and (a, b, c) of the faces */
		get_co_ordinates (vertices, faces, x, y, z, a, b, c, fp);

		if ((vert  = calloc (vertices, sizeof (vec))) == NULL)
		{
				perror ("calloc");
				return ERR;
		}

		for (i = 0; i < vertices; i++)
		{
				vert[i].x = x[i];
				vert[i].y = y[i];
				vert[i].z = z[i];
		}

		/* Find minimum and maximum */
		find_min_max (vertices, x, y, z, &min, &max);

		/* Get center value */
		get_center (min, max, &center);

		/* Get the value of extent */
		e = get_extent (min, max);

		/* Initialize camera and up vectors */
		initialize_vector (&camera, 1, 0, 0);
		initialize_vector (&up, 0, 0, 1);

		/* Get rotation matrix */
		rotation(1, &r_x, theta1); 
		rotation(2, &r_y, theta2); 
		rotation(3, &r_z, theta3); 

		/* Rotate the points */
		matrix_multiplication (r_x, r_y, A, 3);
		matrix_multiplication (A, r_z, B, 3);
		vector_rotation (&camera, &up, B);

		/* Move and scale the camera vector */
		move_scale (e, &camera, &center);		

		/* Determine 3D co-ordinates */
		determine_3D (&left, &right, &up, &bottom, &top, &topleft, &center, camera, e);

		/* Find pixel */
		find_pixel (vert, camera, a, b, c, rows, cols, faces, out_fp, left, right, top, bottom, topleft);

		fcloseall ();

		return SUCCESS;
}

int get_num_ver_faces (int *faces, int *vertices, FILE *fp)
{
		/* Declarations */
		char format[4];
		char data[20];

		fseek (fp, 0, SEEK_SET);

		if (fscanf (fp, "%s", format) <= 0)
		{
				printf ("EOF reached or error\n");
				return ERR;
		}

		if (strcmp (format, "ply") != 0)
		{
				printf ("Input \"ply\" format file\n");
				return ERR;
		}

		while (1)
		{
				if (fscanf (fp, "%s", data) <= 0)
				{
						return ERR;
				}

				if (strcmp (data, "vertex") == 0)
				{
						if (fscanf (fp, "%d", vertices) <= 0)
						{
								return ERR;
						}
						break;
				}
		}

		while (1)
		{
				if (fscanf (fp, "%s", data) <= 0)
				{
						return ERR;
				}

				if (strcmp (data, "face") == 0)
				{
						if (fscanf (fp, "%d", faces) <= 0)
						{
								return ERR;
						}
						break;
				}
		}

		return SUCCESS;
}

int get_co_ordinates (int vertices, int faces, float *x, float *y, float *z, int *a, int *b, int *c, FILE *fp)
{
		int temp, i;
		char data[20];

		fseek (fp, 0, SEEK_SET);

		while (fscanf (fp, "%s", data) > 0)
		{
				if (strcmp (data, "end_header") == 0)
				{
						break;
				}
		}

		for (i = 0; i < vertices; i++)
		{
				if (fscanf (fp, "%f %f %f", &x[i], &y[i], &z[i]) <= 0)
				{
						return ERR;
				}
		}

		for (i = 0; i < faces; i++)
		{
				if (fscanf (fp, "%d %d %d %d", &temp, &a[i], &b[i], &c[i]) <= 0)
				{
						return ERR;
				}
		}

		return SUCCESS;
}

int print_vertices (int vertices, float *x, float *y, float *z)
{
		int i;

		for (i = 0; i < vertices; i++)
		{
				printf ("%f %f %f\n", x[i], y[i], z[i]);
		}
}

int print_faces (int faces, int *a, int *b, int *c)
{
		int i;

		for (i = 0; i < faces; i++)
		{
				printf ("%d %d %d\n", a[i], b[i], c[i]);
		}
}

int find_min_max (int vertices, float *x, float *y, float *z, vec *min, vec *max)
{
		/* Declarations */
		float min_x = x[0], min_y = y[0], min_z = z[0];
		float max_x = x[0], max_y = y[0], max_z = z[0];
		int i;

		for (i = 1; i < vertices; i++)
		{
				if (x[i] < min_x)
				{
						min_x = x[i];
				}

				if (y[i] < min_y)
				{
						min_y = y[i];
				}

				if (z[i] < min_z)
				{
						min_z = z[i];
				}

				if (x[i] > max_x)
				{
						max_x = x[i];
				}

				if (y[i] > max_y)
				{
						max_y = y[i];
				}

				if (z[i] > max_z)
				{
						max_z = z[i];
				}
		}

		min -> x = min_x;
		min -> y = min_y;
		min -> z = min_z;

		max -> x = max_x;
		max -> y = max_y;
		max -> z = max_z;
}

int get_center (vec min, vec max, vec *center)
{
		center -> x = (min.x + max.x)/2;
		center -> y = (min.y + max.y)/2;
		center -> z = (min.z + max.z)/2;
}

float get_extent (vec min, vec max)
{
		float a = max.x - min.x, b = max.y - min.y, c = max.z - min.z;

		if ((a >= b) && (a >= c))
		{
				return a;
		}

		if ((b >= a) && (b >= c))
		{
				return b;
		}


		if ((c >= b) && (c >= a))
		{
				return c;
		}
}

void initialize_vector (vec *camera, int x, int y, int z)
{
		camera->x = (float)x;
		camera->y = (float)y;
		camera->z = (float)z;
}

void rotation(int x, float ***r, float theta)
{
		int i;

		*r = (float **)calloc (3, sizeof (float *));

		for (i = 0; i < 3; i++)
		{
				(*r)[i] = calloc (3, sizeof (int));
		}

		/* Rotation matrix of x-axis */
		if (x == 1)
		{
				(*r)[0][0] = 1;
				(*r)[0][1] = 0;
				(*r)[0][2] = 0;
				(*r)[1][0] = 0;
				(*r)[1][1] = cos (theta);
				(*r)[1][2] = -sin (theta);
				(*r)[2][0] = 0;
				(*r)[2][1] = sin (theta);
				(*r)[2][2] = cos (theta);
		}

		/* Rotation matrix of y-axis */
		else if (x == 2)
		{
				(*r)[0][0] = cos (theta);
				(*r)[0][1] = 0;
				(*r)[0][2] = sin (theta);
				(*r)[1][0] = 0;
				(*r)[1][1] = 1;
				(*r)[1][2] = 0;
				(*r)[2][0] = -sin (theta);
				(*r)[2][1] = 0;
				(*r)[2][2] = cos (theta);
		}

		/* Rotation matrix of z-axis */
		else
		{
				(*r)[0][0] = cos (theta);
				(*r)[0][1] = -sin (theta);
				(*r)[0][2] = 0;
				(*r)[1][0] = sin (theta);
				(*r)[1][1] = cos (theta);
				(*r)[1][2] = 0;
				(*r)[2][0] = 0;
				(*r)[2][1] = 0;
				(*r)[2][2] = 1;
		}
}

void print_matrix (float **r)
{
		int i, j;

		for (i = 0; i < 3; i++)
		{
				for (j = 0; j < 3; j++)
				{
						printf ("%f\t", r[i][j]);
				}
				printf ("\n");
		}
}

float covert_rad (float theta)
{
		return (M_PI / 180.0 * theta);
}

void move_scale (float extent, vec *camera, vec *center)
{
		camera->x = 1.5 * extent * camera->x + center->x;
		camera->y = 1.5 * extent * camera->y + center->y;
		camera->z = 1.5 * extent * camera->z + center->z;
}

void cross_product (vec A, vec B, vec *C)
{
		/* (a * b) = <(ay*bz - az * by), (az * bx - ax * bz), (ax * by - ay * bx)>*/
		C->x = A.y * B.z - A.z * B.y;
		C->y = A.z * B.x - A.x * B.z;
		C->z = A.x * B.y - A.y * B.x;
}

void sum (vec A, vec B, vec *C)
{
		C->x = A.x + B.x;
		C->y = A.y + B.y;
		C->z = A.z + B.z;
}

void sub (vec A, vec B, vec *C)
{
		C->x = A.x - B.x;
		C->y = A.y - B.y;
		C->z = A.z - B.z;
}

void determine_3D (vec *left, vec *right, vec *up, vec *bottom, vec *top, vec *topleft, vec *center, vec camera, float e)
{
		vec temp;
		float a, e_2a, e_2;

		/* <left> = <up> * <center - camera>*/
		sub (*center, camera, &temp);
		cross_product (*up, temp, left);


		/* a = ||<left>|| */
		a = norm (*left);
		e_2a = e / (2.0 * a);
		e_2 = e / 2.0;

		/* <left> = (E/2a) <left> + <center> */
		temp = *left;
		mul (e_2a, &temp);
		sum (temp, *center, left);

		/* <right> = <center - camera> * <up>*/
		sub (*center, camera, &temp);
		cross_product (temp, *up, right);

		/* <right> = E/(2 * a)<right> + <center> */
		temp = *right;
		mul (e_2a, &temp);
		sum (temp, *center, right);

		/* <top> = E/(2)<up> + <center> */
		temp = *up;
		mul (e_2, &temp);
		sum (temp, *center, top);

		/* <bottom> = -(E/2)<up> + <center> */
		temp = *up;
		mul (-e_2, &temp);
		sum (temp, *center, bottom);

		/* <topleft> = E/2 <up> + <left> */
		temp = *up;
		mul (e_2, &temp);
		sum (temp, *left, topleft);

		return;
}

float norm (vec A)
{
		return (sqrt (A.x * A.x + A.y * A.y + A.z * A.z));
}

void mul (float num, vec *A)
{
		A->x = num * A->x;
		A->y = num * A->y;
		A->z = num * A->z;
}


void print_cord (vec left, vec right, vec up, vec bottom, vec topleft, vec center, vec camera)
{
		printf ("left\n");
		print_vector (left);
		printf ("right\n");
		print_vector (right);
		printf (" up\n");
		print_vector (up);
		printf ("bottom\n");
		print_vector (bottom);
		printf ("topleft\n");
		print_vector (topleft);
		printf ("center\n");
		print_vector (center);
		printf ("camera\n");
		print_vector (camera);
}

void print_vector (vec A)
{
		printf ("(x, y, z): (%f, %f, %f)\n", A.x, A.y, A.z);
}

int vector_matrix_multiplication (float *V, float **A, int index, float *R)
{
		/* Declarations */
		int i, j;

		*R = 0;

		for (i = 0; i < 3; i++)
		{
				*R += V[i] * A[i][index]; 
		}

		return SUCCESS;
}

int matrix_multiplication (float **A, float **B, float **C, int num)
{
		/* Declarations */
		int i, j;
		float temp;

		for (i = 0; i < num; i++)
		{
				for (j = 0; j < num; j++)
				{
						vector_matrix_multiplication (A[i], B, j, &temp);
						C[i][j] = temp;
				}
		}
}

void vector_rotation (vec *camera, vec *up, float **A)
{
		int i, j;
		float temp[3], result;


		for (i = 0; i < 2; i++)
		{
				if (i == 0)
				{
						temp[0] = camera -> x;
						temp[1] = camera -> y;
						temp[2] = camera -> z;
				}

				if (i == 1)
				{
						temp[0] = up -> x;
						temp[1] = up -> y;
						temp[2] = up -> z;
				}

				for (j = 0; j < 3; j++)
				{
						vector_matrix_multiplication (temp, A, j, &result);

						if ((j == 0))
						{
								if (i == 0)
										camera -> x = result;
								else if (i == 1)
										up -> x = result;
						}

						if ((j == 1))
						{
								if (i == 0)
										camera -> y = result;
								else if (i == 1)
										up -> y = result;
						}

						if ((j == 2))
						{
								if (i == 0)
										camera -> z = result;
								else if (i == 1)
										up -> z = result;
						}
				}
		}
}

void vector_coordinates (vec *image, vec right, vec left, vec bottom, vec top, vec topleft, int cols, int rows)
{
		int i, j;
		vec temp1, temp2, temp3, temp4;

		sub (right, left, &temp1);
		sub (bottom, top, &temp2);

		mul ((float)cols/(float)255, &temp1);
		mul ((float)rows/(float)255, &temp2);

		sum (topleft, temp1, &temp3);
		sum (temp3, temp2, image);

		return;
}

void dot_product (vec A, vec B, float *C)
{
		int i;

		for (i = 0; i < 3; i++)
		{
				*C = A.x * B.x + A.y * B.y + A.z * B.z;
		}
}

int find_intersect (int vertices, float *n, float *d, int *skip, vec *image, vec camera, vec *intersect)
{
		/* Declarations */
		int i = 0, j = 0;
		vec temp1, temp2;

		for (i = 0; i < vertices; i++)
		{
				if (d[i] < 0.001)
				{
						intersect[i].x = intersect[i].y = intersect[i].z = 0;
						skip[j] = i;
						j++;
						continue;
				}

				sub (image[i], camera, &temp1);
				mul ((n[i] / d[i]), &temp1);
				sum (camera, temp1, (intersect + i));
		}

		return j;
}

void intersection_point (float *dot1, float *dot2, float *dot3, vec intersect, vec v0, vec v1, vec v2)
{
		/* Declarations */
		vec temp1, temp2, temp3, temp4;

		/* Calculate dot1 */
		sub (v2, v0, &temp1);
		sub (v1, v0, &temp2);
		cross_product (temp1, temp2, &temp3);

		sub (intersect, v0, &temp1);
		sub (v1, v0, &temp2);
		cross_product (temp1, temp2, &temp4);

		dot_product (temp3, temp4, dot1);

		/* Calculate dot2 */
		sub (v0, v1, &temp1);
		sub (v2, v1, &temp2);
		cross_product (temp1, temp2, &temp3);

		sub (intersect, v1, &temp1);
		sub (v2, v1, &temp2);
		cross_product (temp1, temp2, &temp4);

		dot_product (temp3, temp4, dot2);

		/* Calculate dot3 */
		sub (v1, v2, &temp1);
		sub (v0, v2, &temp2);
		cross_product (temp1, temp2, &temp3);

		sub (intersect, v2, &temp1);
		sub (v0, v2, &temp2);
		cross_product (temp1, temp2, &temp4);

		dot_product (temp3, temp4, dot3);
}

void find_pixel (vec *vert, vec camera, int *a, int *b, int *c, int rows, int cols, int faces, FILE *fp, vec left, vec right, vec top, vec bottom, vec topleft)
{
		int i, j, k;
		vec temp, temp2, temp3, temp4, temp5, intersect, tmp1, tmp2, tmp3, tmp4, tmp5, ABC, image;
		unsigned char *pixel = calloc (256 * 256, 1);
		float dot1, dot2, dot3, n, d, temp1, D;

		if (fp == NULL)
		{
				printf ("file\n");
				return;
		}
		fseek (fp, 0, SEEK_SET);

		fprintf (fp, "%s %d %d %d ", "P5", 256, 256, 255);

		for (i = 0; i < (rows * cols); i++)
		{
				if ((i % 100) == 0)
						printf ("complete %f\r", (float)i * 100 / (rows * cols));
				/* Calculate image vector */
				vector_coordinates (&image, right, left, bottom, top, topleft, (i % 256), (i / 256));

				for (j = 0; j < faces; j++)
				{
						/* Calulate ABC */
						sub (vert[b[j]], vert[a[j]], &tmp1);
						sub (vert[c[j]], vert[a[j]], &tmp2);
						cross_product (tmp1, tmp2, &ABC);

						/* Calculate D */
						dot_product (ABC, vert[a[j]], &D);
						D = -D;

						/* Calculate n */
						dot_product (ABC, camera, &temp1);
						n = - temp1 - D;

						/* Calculate d */
						sub (image, camera, &temp);
						dot_product (ABC, temp, &d);

						/* If d is near zero, skip */
						if (fabs(d) < 0.0001)
						{
								continue;
						}

						/* Find intersect */
						sub (image, camera, &temp2);
						mul ((n / d), &temp2);
						sum (camera, temp2, &intersect);

						/* Find dot1, dot2 and dot3 */
						intersection_point (&dot1, &dot2, &dot3, intersect, vert[a[j]], vert[b[j]], vert[c[j]]);

						/* Find dot1, dot2 and dot3 */
						if ((dot1 >= 0.0) && (dot2 >= 0.0) && (dot3 >= 0.0))
						{
								pixel[i] = 155 + j % 100;;
						}
				}
		}

		for (i = 0; i < 256 * 256; i++)
		{
				fwrite (&pixel[i], 1, 1, fp);
		}
}
