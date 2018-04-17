#pragma once

#include "glm.h"
#include "dynamic_body.h"
#include "physics_types.h"
#include "physics_math.h"
#include <vector>
#include <minmax.h>

class RigidBody : public DynamicBody
{
public:	
	std::vector<glm::vec3> verts;
	std::vector<PhysicsMath::Triangle> tris;
	
	glm::vec3 rot = glm::vec3(0,0,0);

	float
		mass = 1,
		denisty = 1,
		maxRadius,
		e = 0.5;

	glm::mat3 I;





	RigidBody() : DynamicBody()
	{		
		type = RIGID_BODY;
	}

	void rotate()
	{
		rot.r += rot.b;
	}

	virtual glm::mat3 setInertiaTensor()
	{
		float * ii = new float[9];
		for (int i = 0; i < 9; ++i)
			ii[i] = 0;
		for (int i = 0; i < verts.size(); ++i)
		{
			ii[0] += verts[i].y*verts[i].y + verts[i].z*verts[i].z;
			ii[1] += -(verts[i].x * verts[i].y);
			ii[2] += -(verts[i].x * verts[i].z);
			ii[3] += -(verts[i].x * verts[i].y);
			ii[4] += verts[i].x*verts[i].x + verts[i].z*verts[i].z;
			ii[5] += -(verts[i].y * verts[i].z);
			ii[6] += -(verts[i].x * verts[i].z);
			ii[7] += -(verts[i].y * verts[i].z);
			ii[8] += verts[i].x*verts[i].x + verts[i].y*verts[i].y;
		}
		for (int i = 0; i < 9; ++i)
			ii[i] *= mass / 12;
		return	(I = glm::make_mat3(ii));
	}
	void setMaxR()
	{
		for (int i = 0; i < verts.size(); ++i)
		{
			float r = (pos - verts[i]).length();
			if (r > maxRadius)
				maxRadius = r;
		}
	}

};



class RigidCube : public RigidBody
{
public:
	


};

glm::vec3 checkCollision(RigidBody* b1, RigidBody* b2);
