#pragma once

#include <vector>
#include "dynamic_body.h"
#include "physics_types.h"
#include "rigid_body.h"

struct PhysicsManager
{
	float dt = 1;

	std::vector<DynamicBody *> bodies;

	PhysicsManager(){}
	PhysicsManager(std::vector<DynamicBody *> _bodies) 
	{
		bodies = _bodies;
	}

	DynamicBody * addBody(DynamicBody * body)
	{
		bodies.push_back(body);
		return body;
	}

	void updateBodies()
	{
		for (DynamicBody * b : bodies)
		{
			if(!b->isStatic)
				b->pos += glm::vec3(0,-1,0) * dt;
			b->update(dt);

			if (b->type == RIGID_BODY)
			{
				RigidBody* rb = dynamic_cast<RigidBody*>(b);
				rb->rotate();
			}
		}

		//for (DynamicBody * b1 : bodies)
		//	if (b1->type == RIGID_BODY)
		//		for (DynamicBody * b2 : bodies)
		//			if(b1 != b2)
		//				if (b2->type == RIGID_BODY)
		//				{
		//					RigidBody* rb1 = dynamic_cast<RigidBody*>(b1);
		//					RigidBody* rb2 = dynamic_cast<RigidBody*>(b2);
		//					checkCollision(rb1, rb2);
		//				}

		checkCollision(dynamic_cast<RigidBody*>(bodies[0]), dynamic_cast<RigidBody*>(bodies[1]));

	}

	void setDeltaTime(float _dt)
	{
		dt = _dt;
	}

	RigidBody * getRigidBody(int id)
	{
		if (bodies[id]->type == RIGID_BODY)
			return dynamic_cast<RigidBody*>(bodies[id]);
		else
			return nullptr;
	}

	DynamicBody * getDynamicBody(int id)
	{
		if (bodies[id]->type == DYNAMIC_BODY)
			return bodies[id];
		else
			return nullptr;
	}

	void printOut()
	{
		for (DynamicBody * b : bodies)
			printf("%i\n", b->type);
	}
};
