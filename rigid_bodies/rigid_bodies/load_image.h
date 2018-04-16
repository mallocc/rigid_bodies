#pragma once

#include "opengl.h"

struct image_data
{
	unsigned char * data;
	int w, h, n;
};

GLuint load_texture_from_image(const char *fname);
GLuint load_texture_blank();
GLuint load_texture_uniform(int r, int g, int b);

image_data get_data(const char * filename);