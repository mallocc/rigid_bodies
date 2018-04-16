#pragma once

#include "glm.h"

struct Light
{
	glm::vec3 pos;
	glm::vec3 color;
	float brightness;
	float specular_scale;
	float shininess;
};