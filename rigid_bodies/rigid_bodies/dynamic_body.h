#pragma once

#include "glm.h"
#include "physics_types.h"

class DynamicBody
{
public:
	unsigned int type;
	glm::vec3
		pos,
		vel;

	bool isStatic = 0;

	DynamicBody()
	{
		type = DYNAMIC_BODY;
	}

	virtual void update(float dt) 
	{
		pos += vel * dt;
	};
};