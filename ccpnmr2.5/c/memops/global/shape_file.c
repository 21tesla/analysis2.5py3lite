
/*
======================COPYRIGHT/LICENSE START==========================

shape_file.c: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/
#include "shape_file.h"

Shape_file new_shape_file(int ndim, int ncomponents, int *points)
{
    int i, j, n;
    Shape_file shape_file = (Shape_file) NULL;

    if (ndim > MAX_NDIM)
	return NULL;

    MALLOC_NEW(shape_file, struct Shape_file, 1);
    MALLOC_NEW(shape_file->components, struct Shape_component, ncomponents);

    for (i = 0; i < ncomponents; i++)
    {
	shape_file->components[i].amplitude = 1;
	shape_file->components[i].nshapes = ndim;
	MALLOC_NEW(shape_file->components[i].shapes, struct Shape_shape, ndim);

	for (j = 0; j < ndim; j++)
	{
	    shape_file->components[i].shapes[j].size = points[j];
	    shape_file->components[i].shapes[j].offset = 0;
	    MALLOC_NEW(shape_file->components[i].shapes[j].values, float, points[j]);
	    ZERO_VECTOR(shape_file->components[i].shapes[j].values, points[j]);
	}
    }

    shape_file->ndim = ndim;
    shape_file->ncomponents = ncomponents;
    COPY_VECTOR(shape_file->points, points, ndim);
    CUMULATIVE(shape_file->cum_points, points, n, ndim);

    return shape_file;
}

void delete_shape_file(Shape_file shape_file)
{
    int i, j;

/*
    printf("delete_shape_file\n");
*/

    if (!shape_file)
	return;

    for (i = 0; i < shape_file->ncomponents; i++)
    {
	for (j = 0; j < shape_file->ndim; j++)
	    FREE(shape_file->components[i].shapes[j].values, float);

	FREE(shape_file->components[i].shapes, struct Shape_shape);
    }

    FREE(shape_file->components, struct Shape_component);
    FREE(shape_file, struct Shape_file);
}

CcpnStatus set_amplitude_shape_component(Shape_file shape_file, int comp,
	float amplitude, CcpnString error_msg)
{
    if ((comp < 0) || (comp >= shape_file->ncomponents))
    {
	sprintf(error_msg, "comp = %d, must be in range [0, %d]", comp, shape_file->ncomponents-1);
	return CCPN_ERROR;
    }

    shape_file->components[comp].amplitude = amplitude;

    return CCPN_OK;
}

CcpnStatus set_data_shape_shape(Shape_file shape_file, int comp, int shape,
        int size, int offset, float *values, CcpnString error_msg)
{
    Shape_shape *shape_shape;

    if ((comp < 0) || (comp >= shape_file->ncomponents))
    {
	sprintf(error_msg, "comp = %d, must be in range [0, %d]", comp, shape_file->ncomponents-1);
	return CCPN_ERROR;
    }

    if ((shape < 0) || (shape >= shape_file->ndim))
    {
	sprintf(error_msg, "shape = %d, must be in range [0, %d]", shape, shape_file->ndim-1);
	return CCPN_ERROR;
    }

    if ((size < 1) || (size > shape_file->points[shape]))
    {
	sprintf(error_msg, "size = %d, must be in range [1, %d]", size, shape_file->points[shape]);
	return CCPN_ERROR;
    }

    if ((offset < 0) || (offset > (shape_file->points[shape]-size)))
    {
	sprintf(error_msg, "offset = %d, must be in range [0, %d]", offset, shape_file->points[shape]-size);
	return CCPN_ERROR;
    }

    shape_shape = &(shape_file->components[comp].shapes[shape]);

    ZERO_VECTOR(shape_shape->values, shape_file->points[shape]);
    shape_shape->size = size;
    shape_shape->offset = offset;
    COPY_VECTOR(shape_shape->values+offset, values, size);

    return CCPN_OK;
}
 
static float get_point_shape_component(Shape_component *shape_component, int *point)
{
    int i;
    float value = shape_component->amplitude;
    Shape_shape *shapes = shape_component->shapes;

/*
printf("get_point_shape_component: ampl = %4.3f\n", value);
*/
    for (i = 0; i < shape_component->nshapes; i++)
/*
{
printf("  get_point_shape_component: i = %d: point = %d, value = %4.3f\n", i, point[i], shapes[i].values[point[i]]);
*/
	value *= shapes[i].values[point[i]];
/*
}
printf("get_point_shape_component: value = %4.3f\n", value);
*/

    return value;
}

float get_point_shape_file(Shape_file shape_file, int *point,
				int ncomponents, int *components)
{
    int i;
    float value = 0;
    Shape_component *shape_components = shape_file->components;

    if (components)
    {
	for (i = 0; i < ncomponents; i++)
	    value += get_point_shape_component(&shape_components[components[i]], point);
    }
    else
    {
	for (i = 0; i < shape_file->ncomponents; i++)
	    value += get_point_shape_component(&shape_components[i], point);
    }

    return value;
}

